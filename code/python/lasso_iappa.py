"""Thuc nghiem IAPPA tren bai toan elastic-net / LASSO tong hop.

Ban port Python cua `code/matlab/main_experiment.m`, de nguoi doc khong co
MATLAB van tai lap duoc thi nghiem. Bai toan:

    F(x) = 1/2 ||A x - b||^2 + rho/2 ||x||^2 + mu ||x||_1,      x in R^n.

Toan tu prox_{lambda F} khong co cong thuc dong, nen bai toan con

    Phi_lambda(x; y) = F(x) + 1/(2 lambda) ||x - y||^2

duoc giai xap xi bang FISTA noi. So vong lap noi T_k chinh la "nut van" sinh
sai so co kiem soat: T_k nho -> xap xi tho, T_k lon -> xap xi gan chinh xac.

Tu do sinh ra hai chung nhan sai so hop le cho cung mot diem z:

  * Loai 1: eps1 = sqrt(2 lambda * delta_inner), voi
    delta_inner = Phi_lambda(z; y) - min_x Phi_lambda(x; y);
  * Loai 2: eps2 = eps1 * sqrt((lambda rho + 1) / (lambda rho)), suy ra tu
    tinh loi manh he so rho cua F (Menh de ve moi quan he giua cac xap xi).

Script xuat 3 hinh vao `figures/` va mot bang tong ket vao `results/`.

Chay:  python code/python/lasso_iappa.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "figures"
OUT = ROOT / "results"
FIG.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------------------------
# Tham so, giu dung nhu ban MATLAB
# ----------------------------------------------------------------------------
SEED = 7
N, M, S = 80, 30, 8
MU = 0.10  # he so phat sparsity (L1)
RHO = 0.02  # he so loi manh (ridge)
LAMBDA = 2.0  # tham so prox ngoai

TREF = 200  # so vong lap FISTA noi lam moc "gan dung"
K = 120  # so buoc lap ngoai de ve do thi
TREF_REF = 500  # so lap noi cho moc tham chieu F*, x*
K_REF = 250  # so lap ngoai cho moc tham chieu
NOISE_FLOOR = 1e-11  # san sai so so hoc


def F_obj(A: np.ndarray, b: np.ndarray, x: np.ndarray) -> float:
    """Gia tri ham muc tieu F(x)."""
    return float(
        0.5 * np.sum((A @ x - b) ** 2)
        + 0.5 * RHO * np.sum(x**2)
        + MU * np.sum(np.abs(x))
    )


def soft_threshold(v: np.ndarray, tau: float) -> np.ndarray:
    """Toan tu prox cua chuan 1."""
    return np.sign(v) * np.maximum(np.abs(v) - tau, 0.0)


def inner_fista_trace(
    A: np.ndarray,
    b: np.ndarray,
    y: np.ndarray,
    x0: np.ndarray,
    tmax: int,
    L: float,
) -> tuple[np.ndarray, np.ndarray]:
    """FISTA giai xap xi bai toan con, tra ve toan bo vet lap."""
    z = x0.copy()
    x_km1 = x0.copy()
    t = 1.0

    xtrace = np.zeros((x0.size, tmax))
    objtrace = np.zeros(tmax)

    for j in range(tmax):
        grad = A.T @ (A @ z - b) + RHO * z + (z - y) / LAMBDA
        xk = soft_threshold(z - grad / L, MU / L)

        t_next = (1.0 + np.sqrt(1.0 + 4.0 * t * t)) / 2.0
        z = xk + ((t - 1.0) / t_next) * (xk - x_km1)

        x_km1 = xk
        t = t_next

        xtrace[:, j] = xk
        objtrace[j] = (
            0.5 * np.sum((A @ xk - b) ** 2)
            + 0.5 * RHO * np.sum(xk**2)
            + MU * np.sum(np.abs(xk))
            + np.sum((xk - y) ** 2) / (2.0 * LAMBDA)
        )

    return xtrace, objtrace


def run_outer(
    A: np.ndarray,
    b: np.ndarray,
    x0: np.ndarray,
    k_outer: int,
    tref: int,
    schedule,
    L: float,
) -> dict[str, np.ndarray]:
    """Vong lap ngoai IAPPA, tich luy dong thoi so sach loai 1 va loai 2."""
    x = x0.copy()
    u = x0.copy()
    a_param = 1.0
    eta1 = delta1 = delta2 = 0.0

    conv_factor = np.sqrt((LAMBDA * RHO + 1.0) / (LAMBDA * RHO))

    a_vec = np.zeros(k_outer + 1)
    a_vec[0] = a_param
    eps1v = np.zeros(k_outer)
    eps2v = np.zeros(k_outer)
    delta1v = np.zeros(k_outer + 1)
    delta2v = np.zeros(k_outer + 1)
    fvals = np.zeros(k_outer + 1)
    fvals[0] = F_obj(A, b, x0)
    alphav = np.zeros(k_outer)
    etav = np.zeros(k_outer + 1)
    x_last = x0.copy()

    for k in range(1, k_outer + 1):
        alpha = (
            np.sqrt((a_param * LAMBDA) ** 2 + 4.0 * a_param * LAMBDA)
            - a_param * LAMBDA
        ) / 2.0
        y = (1.0 - alpha) * x + alpha * u

        tk = int(min(max(schedule(k), 1), tref))
        xtrace, objtrace = inner_fista_trace(A, b, y, x, tref, L)

        x_new = xtrace[:, tk - 1]
        d_inner = max(objtrace[tk - 1] - objtrace.min(), 0.0)

        eps1 = np.sqrt(2.0 * LAMBDA * d_inner)
        eps2 = eps1 * conv_factor

        a_new = (1.0 - alpha) * a_param
        u = u - (1.0 / alpha) * (y - x_new)

        eta1 = eta1 + eps1 / alpha
        delta1 = (1.0 - alpha) * delta1 + (alpha * eta1) ** 2 / (2.0 * LAMBDA)
        delta2 = (1.0 - alpha) * delta2 + eps2**2 / (2.0 * LAMBDA)

        x = x_new
        a_param = a_new
        a_vec[k] = a_param
        eps1v[k - 1] = eps1
        eps2v[k - 1] = eps2
        delta1v[k] = delta1
        delta2v[k] = delta2
        alphav[k - 1] = alpha
        etav[k] = eta1
        fvals[k] = F_obj(A, b, x)
        x_last = x

    return {
        "Avec": a_vec,
        "eps1": eps1v,
        "eps2": eps2v,
        "delta1": delta1v,
        "delta2": delta2v,
        "Fvals": fvals,
        "alpha": alphav,
        "eta": etav,
        "x_last": x_last,
    }


def run_classical(
    A: np.ndarray, b: np.ndarray, x0: np.ndarray, k_outer: int, tref: int, L: float
) -> np.ndarray:
    """PPA co dien (khong tang toc) de so sanh."""
    x = x0.copy()
    fvals = np.zeros(k_outer + 1)
    fvals[0] = F_obj(A, b, x0)
    for k in range(1, k_outer + 1):
        xtrace, _ = inner_fista_trace(A, b, x, x, tref, L)
        x = xtrace[:, tref - 1]
        fvals[k] = F_obj(A, b, x)
    return fvals


def run_all() -> dict[str, object]:
    """Chay toan bo thi nghiem, tra ve ket qua + so lan vi pham hai chan.

    Tach rieng khoi main() de test kiem tra duoc cac bat dang thuc ly thuyet
    ma khong can sinh hinh hay ghi file.
    """
    rng = np.random.default_rng(SEED)

    A = rng.standard_normal((M, N)) / np.sqrt(M)
    x_true = np.zeros(N)
    x_true[rng.choice(N, S, replace=False)] = 3.0 * rng.standard_normal(S)
    b = A @ x_true + 0.01 * rng.standard_normal(M)

    x0 = np.zeros(N)
    L = np.linalg.norm(A, 2) ** 2 + RHO + 1.0 / LAMBDA

    print(f"Kich thuoc bai toan: n={N}, m={M}, do thua s={S}")

    # -- moc tham chieu do chinh xac cao cho F*, x* --
    ref = run_outer(A, b, x0, K_REF, TREF_REF, lambda k: TREF_REF, L)
    f_star = ref["Fvals"][-1]
    x_star = ref["x_last"]
    print(f"F* (tham chieu, {K_REF} buoc ngoai x {TREF_REF} lap noi): {f_star:.10e}")

    # -- cac lich trinh lap noi, mo phong cac muc sai so khac nhau --
    schedules = {
        "S1": (lambda k: 1, "$T_k=1$ (cực thô)"),
        "S2": (lambda k: 4, "$T_k=4$ (thô)"),
        "S3": (lambda k: round(2 * np.log2(k + 2)) + 1, r"$T_k \sim \log k$"),
        "S4": (lambda k: round(1.5 * np.sqrt(k + 1)) + 2, r"$T_k \sim \sqrt{k}$"),
        "S5": (lambda k: TREF, "$T_k=T_{\\mathrm{ref}}$ (chính xác)"),
    }

    R: dict[str, dict[str, np.ndarray]] = {}
    for name, (fn, _) in schedules.items():
        R[name] = run_outer(A, b, x0, K, TREF, fn, L)
        print(f"{name}: F(x_K)-F* = {R[name]['Fvals'][-1] - f_star:.4e}")

    f_classical = run_classical(A, b, x0, K, TREF, L)
    print(f"Co dien (khong tang toc): F(x_K)-F* = {f_classical[-1] - f_star:.4e}")

    # -- kiem chung cac chan ly thuyet tren lich trinh S4 --
    a0 = R["S4"]["Avec"][0]
    phi0_gap = (F_obj(A, b, x0) - f_star) + (a0 / 2.0) * np.sum((x_star - x0) ** 2)
    beta_k = R["S4"]["Avec"] / a0
    bound1 = beta_k * phi0_gap + R["S4"]["delta1"]
    bound2 = beta_k * phi0_gap + R["S4"]["delta2"]
    actual = R["S4"]["Fvals"] - f_star

    viol1 = int(np.sum(actual > bound1 + 1e-9))
    viol2 = int(np.sum(actual > bound2 + 1e-9))
    print(f"Vi pham chan loai 1 (S4): {viol1} / {K + 1}")
    print(f"Vi pham chan loai 2 (S4): {viol2} / {K + 1}")

    return {
        "A": A, "b": b, "L": L, "f_star": f_star, "x_star": x_star,
        "R": R, "schedules": schedules, "f_classical": f_classical,
        "actual": actual, "bound1": bound1, "bound2": bound2,
        "beta_k": beta_k, "phi0_gap": phi0_gap,
        "viol1": viol1, "viol2": viol2, "n_steps": K + 1,
    }


def main() -> None:
    res = run_all()
    A, b, f_star = res["A"], res["b"], res["f_star"]
    R, schedules, f_classical = res["R"], res["schedules"], res["f_classical"]
    actual, bound1, bound2 = res["actual"], res["bound1"], res["bound2"]
    viol1, viol2 = res["viol1"], res["viol2"]
    x0 = np.zeros(N)

    clip = lambda v: np.maximum(v, NOISE_FLOOR)  # noqa: E731
    kk = np.arange(1, K + 1)
    k_axis = np.arange(K + 1)

    # -- HINH A: so sanh toc do hoi tu --
    plt.figure(figsize=(8, 6))
    plt.loglog(kk, clip(f_classical[1:] - f_star), "-", color="0.3", lw=1.6,
               label="PPA cổ điển")
    for name, color, style in [
        ("S1", "#a3141c", "-"),
        ("S2", "#d95319", "-."),
        ("S4", "#77ab31", "--"),
        ("S5", "#0072bd", "-"),
    ]:
        plt.loglog(kk, clip(R[name]["Fvals"][1:] - f_star), style, color=color,
                   lw=1.7, label=schedules[name][1])
    f0_gap = F_obj(A, b, x0) - f_star
    plt.loglog(kk, 0.5 * f0_gap / kk, ":", color="k", lw=1.2,
               label=r"tham chiếu $O(1/k)$")
    plt.loglog(kk, 2.0 * f0_gap / kk**2, ":", color="0.5", lw=1.2,
               label=r"tham chiếu $O(1/k^2)$")
    plt.ylim(NOISE_FLOOR / 2, 1e2)
    plt.grid(True, which="both", alpha=0.3)
    plt.xlabel("Bước ngoài $k$")
    plt.ylabel("$F(x_k) - F^*$")
    plt.title("So sánh tốc độ hội tụ")
    plt.legend(loc="lower left", fontsize=9)
    plt.tight_layout()
    plt.savefig(FIG / "lasso-convergence-comparison.png", dpi=200)
    plt.close()

    # -- HINH B: chan ly thuyet vs sai so thuc --
    plt.figure(figsize=(8, 6))
    plt.loglog(k_axis[1:], clip(actual[1:]), "-", color="#77ab31", lw=1.8,
               label=r"$F(x_k)-F^*$ thực tế ($T_k\sim\sqrt{k}$)")
    plt.loglog(k_axis[1:], clip(bound1[1:]), "--", color="#d95319", lw=1.5,
               label="chặn từ phân tích loại 1")
    plt.loglog(k_axis[1:], clip(bound2[1:]), "-.", color="#0072bd", lw=1.5,
               label="chặn từ phân tích loại 2")
    plt.grid(True, which="both", alpha=0.3)
    plt.xlabel("Bước ngoài $k$")
    plt.ylabel("giá trị")
    plt.title("Kiểm chứng số học các chặn hội tụ")
    plt.legend(loc="lower left", fontsize=9)
    plt.tight_layout()
    plt.savefig(FIG / "lasso-bound-verification.png", dpi=200)
    plt.close()

    # -- HINH C: suy giam cua eps_k theo tung lich trinh --
    plt.figure(figsize=(8, 6))
    for name, color in [
        ("S1", "#a3141c"),
        ("S2", "#d95319"),
        ("S3", "#edb120"),
        ("S4", "#77ab31"),
        ("S5", "#0072bd"),
    ]:
        plt.loglog(kk, clip(R[name]["eps1"]), "-", color=color, lw=1.4,
                   label=schedules[name][1])
    plt.grid(True, which="both", alpha=0.3)
    plt.xlabel("Bước ngoài $k$")
    plt.ylabel(r"$\varepsilon_k$ loại 1 (sàn tại $10^{-11}$)")
    plt.title("Độ chính xác loại 1 đạt được theo từng lịch trình lặp nội")
    plt.legend(loc="lower left", fontsize=9)
    plt.tight_layout()
    plt.savefig(FIG / "lasso-type1-accuracy.png", dpi=200)
    plt.close()

    summary = pd.DataFrame(
        {
            "schedule": ["PPA cổ điển"] + list(schedules),
            "description": ["không tăng tốc, T=T_ref"]
            + [
                "T_k=1", "T_k=4", "T_k ~ log k", "T_k ~ sqrt(k)", "T_k = T_ref",
            ],
            "final_gap": [f_classical[-1] - f_star]
            + [R[n]["Fvals"][-1] - f_star for n in schedules],
            "final_eps1": [np.nan] + [R[n]["eps1"][-1] for n in schedules],
        }
    )
    summary.to_csv(OUT / "lasso-summary.csv", index=False)

    pd.DataFrame(
        {
            "k": k_axis,
            "F_gap_S4": actual,
            "bound_type1": bound1,
            "bound_type2": bound2,
        }
    ).to_csv(OUT / "lasso-bounds.csv", index=False)

    print()
    print(summary.to_string(index=False))
    print(f"\nF* = {f_star:.10e}   (vi pham chan: loai 1 = {viol1}, loai 2 = {viol2})")


if __name__ == "__main__":
    main()
