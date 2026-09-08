from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "figures"
OUT = ROOT / "results"
FIG.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)

N = 2000
LAM = 1.0
A0 = 1.0
X0 = 4.0
C = 0.20
Q = 1.75


def F(x: float | np.ndarray) -> float | np.ndarray:
    return 0.25 * x**4


def F_conj(v: float | np.ndarray) -> float | np.ndarray:
    return 0.75 * np.abs(v) ** (4.0 / 3.0)


def prox_quartic(y: float, lam: float = 1.0) -> float:
    """Solve z^3 + (z-y)/lam = 0 by Newton's method."""
    z = float(np.cbrt(y)) if abs(y) > 1 else y / (1 + lam)
    for _ in range(80):
        g = z**3 + (z - y) / lam
        gp = 3 * z * z + 1 / lam
        z_new = z - g / gp
        if abs(z_new - z) <= 1e-14 * (1 + abs(z)):
            return float(z_new)
        z = z_new
    return float(z)


def phi(z: float, y: float, lam: float) -> float:
    return float(F(z) + (z - y) ** 2 / (2 * lam))


def gap_type2(z: float, y: float, lam: float) -> float:
    v = (y - z) / lam
    return float(F(z) + F_conj(v) - v * z)


def boundary_points(y: float, lam: float, eps: float, kind: int) -> tuple[float, float]:
    """Find the two boundary points for the type-1 or type-2 certificate."""
    p = prox_quartic(y, lam)
    target = eps**2 / (2 * lam)
    if target == 0:
        return p, p

    points: list[float] = []
    for direction in (-1.0, 1.0):
        def residual(t: float) -> float:
            z = p + direction * t
            value = (
                phi(z, y, lam) - phi(p, y, lam)
                if kind == 1
                else gap_type2(z, y, lam)
            )
            return value - target

        lo, hi = 0.0, max(1e-8, 2 * eps)
        while residual(hi) < 0:
            hi *= 2
        for _ in range(90):
            mid = (lo + hi) / 2
            if residual(mid) >= 0:
                hi = mid
            else:
                lo = mid
        points.append(p + direction * hi)
    return points[0], points[1]


def choose_adversarial(points: tuple[float, float]) -> float:
    """Choose the boundary point having the larger objective value."""
    return float(points[int(np.argmax([F(z) for z in points]))])


def alpha_from_A(A: float, lam: float) -> float:
    c = A * lam
    return (-c + math.sqrt(c * c + 4 * c)) / 2


def run_exact() -> np.ndarray:
    x = nu = X0
    A = A0
    values = [F(x)]
    for _ in range(N):
        a = alpha_from_A(A, LAM)
        y = (1 - a) * x + a * nu
        z = prox_quartic(y, LAM)
        nu = nu - a / ((1 - a) * A * LAM) * (y - z)
        A = (1 - a) * A
        x = z
        values.append(F(x))
    return np.asarray(values, dtype=float)


def run_iappa1(q: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = u = X0
    A = A0
    eta = delta = 0.0
    values = [F(x)]
    cert: list[float] = []
    deltas = [0.0]
    for k in range(N):
        eps = C / (k + 1) ** q
        a = alpha_from_A(A, LAM)
        y = (1 - a) * x + a * u
        z = choose_adversarial(boundary_points(y, LAM, eps, 1))
        p = prox_quartic(y, LAM)
        target = eps**2 / (2 * LAM)
        cert.append((phi(z, y, LAM) - phi(p, y, LAM)) / target)
        u = u - (y - z) / a
        eta = eta + eps / a
        delta = (1 - a) * delta + (a * eta) ** 2 / (2 * LAM)
        A = (1 - a) * A
        x = z
        values.append(F(x))
        deltas.append(delta)
    return np.asarray(values), np.asarray(cert), np.asarray(deltas)


def run_iappa2(q: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = nu = X0
    A = A0
    delta = 0.0
    values = [F(x)]
    cert: list[float] = []
    deltas = [0.0]
    for k in range(N):
        eps = C / (k + 1) ** q
        a = alpha_from_A(A, LAM)
        y = (1 - a) * x + a * nu
        z = choose_adversarial(boundary_points(y, LAM, eps, 2))
        target = eps**2 / (2 * LAM)
        cert.append(gap_type2(z, y, LAM) / target)
        nu = nu - a / ((1 - a) * A * LAM) * (y - z)
        delta = (1 - a) * delta + eps**2 / (2 * LAM)
        A = (1 - a) * A
        x = z
        values.append(F(x))
        deltas.append(delta)
    return np.asarray(values), np.asarray(cert), np.asarray(deltas)


def log_slope(values: np.ndarray, start: int = 500, end: int = 2000) -> float:
    idx = np.arange(start, end + 1)
    mask = np.isfinite(values[idx]) & (values[idx] > 0)
    return float(np.polyfit(np.log(idx[mask]), np.log(values[idx][mask]), 1)[0])


def main() -> None:
    exact = run_exact()
    iappa1, cert1, delta1 = run_iappa1(Q)
    iappa2, cert2, delta2 = run_iappa2(Q)

    k = np.arange(N + 1)
    pd.DataFrame(
        {
            "k": k,
            "exact": exact,
            "IAPPA1_q1.75": iappa1,
            "IAPPA2_q1.75": iappa2,
            "delta1_bound_component": delta1,
            "delta2_bound_component": delta2,
        }
    ).to_csv(OUT / "quartic-results.csv", index=False)

    slope1 = log_slope(delta1)
    slope2 = log_slope(delta2)
    dev1 = float(np.max(np.abs(cert1 - 1.0)))
    dev2 = float(np.max(np.abs(cert2 - 1.0)))

    summary = pd.DataFrame(
        {
            "method": ["Exact accelerated", "IAPPA1 q=7/4", "IAPPA2 q=7/4"],
            "final_objective": [exact[-1], iappa1[-1], iappa2[-1]],
            "delta_log_slope_500_2000": [np.nan, slope1, slope2],
            "max_abs_certificate_minus_1": [np.nan, dev1, dev2],
        }
    )
    summary.to_csv(OUT / "quartic-summary.csv", index=False)

    plt.figure(figsize=(8, 5))
    plt.loglog(k[1:], exact[1:], label="Thuật toán tham chiếu")
    plt.loglog(k[1:], iappa1[1:], label="IAPPA1 (loại 1)")
    plt.loglog(k[1:], iappa2[1:], label="IAPPA2 (loại 2)")
    plt.xlabel("Số vòng lặp k")
    plt.ylabel("Sai số giá trị hàm F(x_k) - F*")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "quartic-objective-convergence.png", dpi=200)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.loglog(k[1:], delta1[1:], label="IAPPA1: thành phần sai số")
    plt.loglog(k[1:], delta2[1:], label="IAPPA2: thành phần sai số")
    plt.loglog(k[10:], delta1[10] * (k[10:] / 10) ** (-0.5), "--", label="Tham chiếu k^{-1/2}")
    plt.loglog(k[10:], delta2[10] * (k[10:] / 10) ** (-2.0), "--", label="Tham chiếu k^{-2}")
    plt.xlabel("Số vòng lặp k")
    plt.ylabel("Thành phần sai số tích lũy")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "quartic-error-accumulation.png", dpi=200)
    plt.close()

    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
