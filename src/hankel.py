import torch


def dst4(x: torch.Tensor) -> torch.Tensor:
    """DST-IV não normalizada para um tensor 1D."""
    n = x.shape[0]

    y = torch.stack(
        (torch.zeros_like(x), x),
        dim=1,
    ).flatten()

    extension = torch.cat(
        (
            y,
            x.new_zeros(4 * n + 1),
            -torch.flip(y[1:], dims=(0,)),
        )
    )

    return -torch.fft.rfft(extension).imag[1 : 2 * n : 2]


def idst4(x: torch.Tensor) -> torch.Tensor:
    """Inversa matemática da DST-IV não normalizada."""
    return dst4(x) / (2 * x.shape[0])


def hankel_forward(
    fr: torch.Tensor,
    r: torch.Tensor,
    k: torch.Tensor,
    dr: float,
) -> torch.Tensor:
    return 2.0 * torch.pi * dr * dst4(fr * r) / k


def hankel_inverse(
    fk: torch.Tensor,
    r: torch.Tensor,
    k: torch.Tensor,
    dk: float,
) -> torch.Tensor:
    return dk / (4.0 * torch.pi**2) * dst4(fk * k) / r