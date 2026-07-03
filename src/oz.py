import numpy as np
from hankel import hankel_forward, hankel_inverse
from optimizer import AndersonMix, PicardMix
from typing import Tuple
import torch

def uhs(r,params):
    sigma = params
    return torch.where(r <= sigma, torch.tensor(float('inf')), torch.tensor(0.0))


class PyOZ:
    def __init__(self, u: callable, params: np.ndarray, closure: str = 'HNC', device = 'cpu', dtype = torch.float64):
        """
        u      : pair potential function u(r, params)
        params : parameters for u
        """
        self.u = u
        self.params = params
        
        # closure: 'PY', 'HNC', 'MSA'
        self.closure = closure
        # device for torch tensors
        self.device = device
        self.dtype = dtype


    def set_closure(self, closure: str = 'PY'):
        if closure not in ('PY', 'HNC', 'MSA'):
            raise ValueError("Closure must be 'PY', 'HNC', or 'MSA'")
        self.closure = closure

    def _closure_relation(self, gamma: torch.Tensor, r: torch.Tensor) -> torch.Tensor:
        """Compute c(r) from gamma(r) using the selected closure."""
        beta_u = self.u(r, self.params) / self.kBT
        if self.closure == 'PY':
            b_r = torch.log1p(gamma) - gamma
            c_r = torch.exp(-beta_u + gamma + b_r) - gamma - 1
        elif self.closure == 'HNC':
            c_r = torch.exp(-beta_u + gamma) - gamma - 1
        elif self.closure == 'MSA':
            b_r = torch.log1p(gamma) - gamma
            c_r = torch.exp(-beta_u + gamma + b_r) - gamma - 1
            c_r[beta_u < float('inf')] = -beta_u[beta_u < float('inf')]  # c(r) = -beta*u(r) for r > sigma
        else:  # not implemented
            raise NotImplementedError(f"Closure '{self.closure}' is not implemented.")
        return c_r

    def solve(self, rho: float, kBT: float = 1.0, rmax: float = 10.0, dr: float = 0.01, max_iter: int = 1000, method: str = 'Anderson', mix_depth: int = 5, alpha: float = 0.6, atol: float = 1e-6, logoutput: bool = False) -> Tuple[np.ndarray, np.ndarray]:
        """
        Solve the Ornstein-Zernike equation for a given radial grid r.
        Returns (r, h_r, c_r).
        """
        self.rho = rho
        self.kBT = kBT

        Npts = int(rmax / dr)
        dk = np.pi / (Npts * dr)
        i = torch.arange(Npts, device=self.device, dtype=self.dtype)
        r = (i + 0.5) * dr  # avoid r=0 to prevent singularity
        k = (i + 0.5) * dk

        # initial outer functions
        c_r = torch.zeros_like(r, device=self.device)
        h_r = torch.zeros_like(r, device=self.device)
        gamma_r = h_r - c_r  # indirect correlation function

        # Anderson mixer
        if method == 'Anderson':
            mixer = AndersonMix(m=mix_depth, alpha=alpha)
        else:
            mixer = PicardMix(alpha=alpha)

        for it in range(1, max_iter + 1):
            # closure
            c_r[:] = self._closure_relation(gamma_r, r)
            # OZ in k-space
            c_k = hankel_forward(c_r, r, k, dr)
            # test if rho*ck is smaller than 1 to avoid divergence using clamp
            denominator = 1.0 / (1.0 - self.rho * c_k)
            denominator = torch.clamp(denominator, min=0.0, max=1e9)  # avoid divergence
            gamma_k = self.rho * c_k * c_k  * denominator
            # back to r-space
            gamma_new = hankel_inverse(gamma_k, r, k, dk)  # inverse DST-I normalization
             # avoid negative values in g(r)
            c_r[:] = self._closure_relation(gamma_new, r)
            h_r[:] = gamma_new + c_r[:]
            gamma_new[:] = torch.clamp(h_r, min=-1.0) - c_r

            # testing convergence
            residue = (gamma_new - gamma_r)*r
            err = torch.sqrt(dr*torch.pow(residue, 2.0).sum()) / atol

            if logoutput and it % 100 == 0:
                print(f"Iter {it}: error={err.item():.2e}")
            if err < 1.0 and it > mix_depth:
                break

            # mixing
            gamma_r[:] = mixer.update(gamma_new, gamma_r)

        self.Niter = it # store the number of iterations
        
        c_r[:] = self._closure_relation(gamma_r, r)
        h_r[:] = gamma_r + c_r
        return r, h_r, c_r