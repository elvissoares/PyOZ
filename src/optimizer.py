from collections import deque
import torch

class PicardMix:
    """
    Picard mixing for fixed-point iterations:

        x_new = G(x_old)
    """
    
    def __init__(self, alpha: float = 0.5):
        if not 0.0 < alpha <= 1.0:
            raise ValueError("alpha must be in the range (0, 1].")
        self.alpha = alpha

    @torch.no_grad()
    def update(self, x_new: torch.Tensor, x_old: torch.Tensor) -> torch.Tensor:
        """
        Compute the next estimate using Picard mixing.

        Parameters
        ----------
        x_new
            Result of applying the fixed-point map G(x_old).

        x_old
            Estimate from the previous iteration.

        Returns
        -------
        torch.Tensor
            Next estimate, on the same device, dtype, and shape as x_old.
        """
        return x_old + self.alpha * (x_new - x_old)


class AndersonMix:
    """
    Anderson mixing para iterações de ponto fixo:

        x_new = G(x_old)

    Os coeficientes c_i são determinados minimizando:

        || sum_i c_i F_i ||²

    sujeito a:

        sum_i c_i = 1

    onde:

        F_i = x_new_i - x_old_i
    """

    def __init__(
        self,
        m: int = 5,
        alpha: float = 1.0e-3,
        regularization: float = 1.0e-8,
    ):
        if m < 1:
            raise ValueError("m deve ser maior ou igual a 1.")

        if not 0.0 < alpha <= 1.0:
            raise ValueError("alpha deve estar no intervalo (0, 1].")

        if regularization < 0.0:
            raise ValueError("regularization deve ser não negativa.")

        self.m = int(m)
        self.alpha = alpha
        self.regularization = regularization

        self.F_hist: deque[torch.Tensor] = deque(maxlen=self.m)
        self.X_hist: deque[torch.Tensor] = deque(maxlen=self.m)

    def reset(self) -> None:
        """Remove todo o histórico do otimizador."""
        self.F_hist.clear()
        self.X_hist.clear()

    @torch.no_grad()
    def update(
        self,
        x_new: torch.Tensor,
        x_old: torch.Tensor,
    ) -> torch.Tensor:
        """
        Calcula a próxima estimativa usando Anderson mixing.

        Parameters
        ----------
        x_new
            Resultado da aplicação do mapa de ponto fixo G(x_old).

        x_old
            Estimativa da iteração anterior.

        Returns
        -------
        torch.Tensor
            Próxima estimativa, no mesmo device, dtype e formato de x_old.
        """

        self._validate_inputs(x_new, x_old)

        residual = x_new - x_old

        # detach evita armazenar o grafo computacional de todas as iterações.
        self.F_hist.append(residual.detach().clone())
        self.X_hist.append(x_new.detach().clone())

        # Mistura linear enquanto ainda não existe histórico suficiente.
        if len(self.F_hist) < self.m:
            return x_old + self.alpha * residual

        # Cada residual pode ser multidimensional. Para construir F, cada
        # residual é convertido para um vetor de tamanho N.
        #
        # F possui formato (N, m).
        F = torch.stack(
            [f.reshape(-1) for f in self.F_hist],
            dim=1,
        )

        # torch.linalg.solve não é normalmente implementado para float16
        # ou bfloat16. A álgebra do pequeno sistema m x m é feita em
        # float32, preservando float64 quando esse for o dtype original.
        solve_dtype = (
            torch.float64
            if F.dtype == torch.float64
            else torch.float32
        )

        F_solve = F.to(dtype=solve_dtype)

        # Matriz de Gram: C_ij = <F_i, F_j>
        C = F_solve.T @ F_solve

        history_size = C.shape[0]

        identity = torch.eye(
            history_size,
            device=C.device,
            dtype=C.dtype,
        )

        # Regularização relativa à magnitude dos resíduos.
        diagonal_scale = torch.diagonal(C).mean()
        machine_epsilon = torch.finfo(C.dtype).eps

        diagonal_scale = torch.clamp(
            diagonal_scale,
            min=machine_epsilon,
        )

        C_regularized = (
            C
            + self.regularization * diagonal_scale * identity
        )

        ones = torch.ones(
            history_size,
            device=C.device,
            dtype=C.dtype,
        )

        try:
            # Resolve:
            #
            #     C u = 1
            #
            # seguido por:
            #
            #     c = u / sum(u)
            #
            # Isso impõe sum(c_i) = 1.
            coefficients = torch.linalg.solve(
                C_regularized,
                ones,
            )

            denominator = coefficients.sum()

            if (
                not torch.isfinite(coefficients).all()
                or not torch.isfinite(denominator)
                or denominator.abs() <= machine_epsilon
            ):
                return x_old + self.alpha * residual

            coefficients = coefficients / denominator

        except RuntimeError:
            # Inclui matrizes singulares ou falhas do backend de álgebra
            # linear.
            return x_old + self.alpha * residual

        # Retorna os coeficientes para o dtype original.
        coefficients = coefficients.to(dtype=x_new.dtype)

        # X possui formato:
        #
        #     (m, *x_new.shape)
        X = torch.stack(list(self.X_hist), dim=0)

        # Converte coefficients de (m,) para:
        #
        #     (m, 1, 1, ...)
        #
        # permitindo broadcasting para qualquer dimensão de x_new.
        coefficient_shape = (
            history_size,
            *([1] * x_new.ndim),
        )

        coefficients = coefficients.reshape(coefficient_shape)

        return torch.sum(coefficients * X, dim=0)

    @staticmethod
    def _validate_inputs(
        x_new: torch.Tensor,
        x_old: torch.Tensor,
    ) -> None:
        if not isinstance(x_new, torch.Tensor):
            raise TypeError("x_new deve ser um torch.Tensor.")

        if not isinstance(x_old, torch.Tensor):
            raise TypeError("x_old deve ser um torch.Tensor.")

        if x_new.shape != x_old.shape:
            raise ValueError(
                "x_new e x_old devem possuir o mesmo formato. "
                f"Recebidos: {x_new.shape} e {x_old.shape}."
            )

        if x_new.device != x_old.device:
            raise ValueError(
                "x_new e x_old devem estar no mesmo device. "
                f"Recebidos: {x_new.device} e {x_old.device}."
            )

        if x_new.dtype != x_old.dtype:
            raise ValueError(
                "x_new e x_old devem possuir o mesmo dtype. "
                f"Recebidos: {x_new.dtype} e {x_old.dtype}."
            )

        if not x_new.is_floating_point():
            raise TypeError(
                "x_new e x_old devem possuir dtype de ponto flutuante."
            )