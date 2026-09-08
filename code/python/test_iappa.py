"""Kiem chung so hoc cho hai thuc nghiem IAPPA.

Muc dich cua tep nay KHONG phai la kiem tra "script co chay khong" -- viec do
`make figures` da lam roi. Muc dich la kiem tra **ket luan toan hoc van con
dung**: neu ai do sua nham mot dau trong cong thuc cap nhat, hoac lam hong
cach dung diem xap xi, cac phep khang dinh duoi day phai do.

Moi phep khang dinh deu ung voi mot phat bieu cu the trong bao cao:

  * `prox` va `soft_threshold` giai dung dieu kien toi uu bac mot;
  * `F_conj` dung la lien hop Fenchel cua F (bat dang thuc Fenchel-Young,
    dat dau bang tai v = F'(x));
  * cac diem duoc dung nam DUNG tren bien ngan sach sai so (chung chi
    rho_k ~ 1), tuc khong vo tinh tao mot nhieu co loi;
  * bac triet tieu do duoc cua delta_k khop bac ly thuyet: -1/2 cho IAPPA1
    va -2 cho IAPPA2 (Dinh ly 3.5 / 3.8 cua bao cao);
  * hai chan hoi tu KHONG bi vi pham o bat ky buoc nao;
  * he so chuyen chung nhan loai 1 -> loai 2 dung bang
    sqrt((lambda*rho + 1)/(lambda*rho)).

Chay:  pytest code/python -q
"""

from __future__ import annotations

import numpy as np
import pytest

import lasso_iappa as lasso
import quartic_iappa as quartic

# ---------------------------------------------------------------------------
# Chay moi thuc nghiem dung MOT lan cho ca module (lasso mat ~10s)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def quartic_res():
    return quartic.run_all()


@pytest.fixture(scope="module")
def lasso_res():
    return lasso.run_all()


# ===========================================================================
# 1. Cac khoi xay dung co ban
# ===========================================================================


@pytest.mark.parametrize("y", [-8.0, -1.0, -0.25, 0.0, 0.25, 1.0, 8.0, 40.0])
@pytest.mark.parametrize("lam", [0.5, 1.0, 3.0])
def test_prox_quartic_thoa_dieu_kien_toi_uu(y: float, lam: float) -> None:
    """z = prox_{lam F}(y) voi F(x)=x^4/4  <=>  z^3 + (z-y)/lam = 0."""
    z = quartic.prox_quartic(y, lam)
    residual = z**3 + (z - y) / lam
    assert abs(residual) < 1e-9 * max(1.0, abs(y))


@pytest.mark.parametrize("y", [-3.0, -0.5, 0.0, 0.5, 3.0])
def test_prox_quartic_la_diem_cuc_tieu_that(y: float) -> None:
    """Phi(z) phai nho hon Phi tai moi diem lan can."""
    lam = 1.0
    z = quartic.prox_quartic(y, lam)
    phi_z = quartic.phi(z, y, lam)
    for h in (-1e-3, -1e-4, 1e-4, 1e-3):
        assert quartic.phi(z + h, y, lam) >= phi_z - 1e-15


@pytest.mark.parametrize("x", [-2.0, -0.3, 0.3, 2.0])
def test_fenchel_young_dat_dau_bang_tai_gradient(x: float) -> None:
    """F(x) + F*(v) >= v*x, dat dau bang khi v = F'(x) = x^3."""
    v_star = x**3
    gap = quartic.F(x) + quartic.F_conj(v_star) - v_star * x
    assert abs(gap) < 1e-12

    for v in (v_star - 0.7, v_star + 0.7):
        assert quartic.F(x) + quartic.F_conj(v) - v * x > -1e-12


@pytest.mark.parametrize("tau", [0.05, 0.4, 1.5])
def test_soft_threshold_la_prox_cua_chuan_1(tau: float) -> None:
    """x = prox_{tau||.||_1}(v)  <=>  (v-x)/tau thuoc duoi vi phan cua |.| tai x."""
    rng = np.random.default_rng(0)
    v = rng.standard_normal(200) * 2.0
    x = lasso.soft_threshold(v, tau)
    sub = (v - x) / tau
    nz, z = np.abs(x) > 1e-12, np.abs(x) <= 1e-12
    # tai toa do khac 0: duoi vi phan la dau cua x
    assert np.allclose(sub[nz], np.sign(x[nz]), atol=1e-12)
    # tai toa do bang 0: duoi vi phan la [-1, 1]
    assert np.all(np.abs(sub[z]) <= 1.0 + 1e-12)


# ===========================================================================
# 2. Ham bac bon: chung chi sai so va bac hoi tu
# ===========================================================================


def test_diem_xap_xi_nam_dung_tren_bien_ngan_sach(quartic_res) -> None:
    """rho_k = (do lech do duoc)/(ngan sach eps_k^2/2lambda) phai bang 1.

    Neu rho_k < 1 thi thuc nghiem da tu cho minh mot xap xi tot hon muc sai so
    khai bao -- ket qua se dep len mot cach khong trung thuc.
    """
    assert quartic_res["max_cert_dev1"] < 1e-6
    assert quartic_res["max_cert_dev2"] < 1e-6
    # va khong duoc lech ve phia "de hon"
    assert np.all(quartic_res["cert1"] <= 1.0 + 1e-6)
    assert np.all(quartic_res["cert2"] <= 1.0 + 1e-6)


def test_bac_triet_tieu_cua_delta_khop_ly_thuyet(quartic_res) -> None:
    """Voi q = 7/4: delta_k ~ k^{-1/2} cho IAPPA1 va k^{-2} cho IAPPA2.

    Day la khac biet trung tam giua sai so loai 1 va loai 2.
    """
    assert quartic_res["slope1"] == pytest.approx(-0.5, abs=0.10)
    assert quartic_res["slope2"] == pytest.approx(-2.0, abs=0.10)
    # hai bac phai tach nhau ro rang, khong chi la sai so do dac
    assert quartic_res["slope1"] - quartic_res["slope2"] > 1.0


def test_ca_ba_phuong_phap_deu_giam_ve_gia_tri_toi_uu(quartic_res) -> None:
    """F_* = 0 nen gia tri cuoi phai rat nho, va nho hon nhieu gia tri dau."""
    for key in ("exact", "iappa1", "iappa2"):
        seq = quartic_res[key]
        assert seq[0] == pytest.approx(quartic.F(quartic.X0))
        assert seq[-1] >= 0.0
        assert seq[-1] < 1e-8 * seq[0]


def test_iappa2_ket_thuc_chinh_xac_hon_iappa1(quartic_res) -> None:
    """Voi cung mot lich trinh sai so, loai 2 phai ve dich sat hon han."""
    assert quartic_res["iappa2"][-1] < quartic_res["iappa1"][-1]


# ===========================================================================
# 3. LASSO: hai chan ly thuyet
# ===========================================================================


def test_hai_chan_hoi_tu_khong_bi_vi_pham(lasso_res) -> None:
    """Khang dinh cot loi cua bao cao.

    Chan cua Dinh ly 2.2:  F(x_k) - F* <= beta_k (phi_0(x*) - F*) + delta_k,
    voi delta_k lay theo so sach loai 1 va loai 2. Truoc day script chi IN ra
    "0/121 vi pham" -- nay neu con dung mot buoc vi pham thi test do.
    """
    assert lasso_res["viol1"] == 0, "chan loai 1 bi vi pham"
    assert lasso_res["viol2"] == 0, "chan loai 2 bi vi pham"
    # kiem tra truc tiep, khong tin vao bien dem
    assert np.all(lasso_res["actual"] <= lasso_res["bound1"] + 1e-9)
    assert np.all(lasso_res["actual"] <= lasso_res["bound2"] + 1e-9)


def test_delta2_khop_cong_thuc_dong_cua_bao_cao(lasso_res) -> None:
    """delta_k = (beta_k/2) * sum_i eps_i^2 / (lambda_i * beta_{i+1}).

    Day la nghiem tuong minh cua truy hoi
    delta_{k+1} = (1-alpha_k) delta_k + eps_k^2/(2 lambda_k)
    duoc phat bieu trong essay/sections/s4 (eq:iappa2-error). Phep kiem nay
    SAC hon test chan hoi tu: chan con nhieu du dia nen mot he so sai van co
    the khong lam no bi vi pham, con dang thuc nay thi sai la do ngay.
    """
    for name in ("S1", "S3", "S4", "S5"):
        r = lasso_res["R"][name]
        beta = r["Avec"] / r["Avec"][0]
        eps = r["eps2"]
        terms = eps**2 / (lasso.LAMBDA * beta[1:])
        expected = beta[1:] / 2.0 * np.cumsum(terms)
        assert np.allclose(r["delta2"][1:], expected, rtol=1e-10, atol=1e-18), name


def test_delta1_khop_cong_thuc_dong_cua_bao_cao(lasso_res) -> None:
    """delta_k = (A beta_k / 2) * sum_i eta_{i+1}^2  (eq:iappa1-error).

    Dang nay chi dung khi quy tac tham so alpha_i^2 = (1-alpha_i) A_i lambda_i
    duoc ton trong -- nen no kiem luon ca quy tac cap nhat alpha.
    """
    for name in ("S1", "S3", "S4", "S5"):
        r = lasso_res["R"][name]
        a0 = r["Avec"][0]
        beta = r["Avec"] / a0
        expected = a0 * beta[1:] / 2.0 * np.cumsum(r["eta"][1:] ** 2)
        assert np.allclose(r["delta1"][1:], expected, rtol=1e-10, atol=1e-18), name


def test_quy_tac_tham_so_alpha_duoc_ton_trong(lasso_res) -> None:
    """alpha_k^2 = (1 - alpha_k) A_k lambda_k, tuc ti so bang dung 1."""
    for name in ("S1", "S4", "S5"):
        r = lasso_res["R"][name]
        alpha, avec = r["alpha"], r["Avec"]
        ratio = alpha**2 / ((1.0 - alpha) * avec[:-1] * lasso.LAMBDA)
        assert np.allclose(ratio, 1.0, rtol=1e-10), name
        assert np.all((alpha >= 0.0) & (alpha < 1.0))


def test_eta_la_tong_tich_luy_cua_eps_chia_alpha(lasso_res) -> None:
    """eta_k = sum_{i<k} eps_i / alpha_i."""
    for name in ("S1", "S3", "S5"):
        r = lasso_res["R"][name]
        expected = np.cumsum(r["eps1"] / r["alpha"])
        assert np.allclose(r["eta"][1:], expected, rtol=1e-12), name


def test_chan_khong_tam_thuong(lasso_res) -> None:
    """Chan phai huu han va that su chan tu tren, khong phai vo cung."""
    assert np.all(np.isfinite(lasso_res["bound1"]))
    assert np.all(np.isfinite(lasso_res["bound2"]))
    assert lasso_res["bound1"][-1] < 1e3
    assert lasso_res["bound2"][-1] < 1e3


def test_beta_k_triet_tieu_bac_hai(lasso_res) -> None:
    """Voi lambda_k hang so, beta_k = O(1/k^2) theo Bo de 2.5."""
    beta = lasso_res["beta_k"]
    assert beta[0] == pytest.approx(1.0)
    assert np.all(np.diff(beta) <= 1e-15)  # khong tang
    k = np.arange(1, len(beta))
    slope = np.polyfit(np.log(k[10:]), np.log(beta[1:][10:]), 1)[0]
    assert slope == pytest.approx(-2.0, abs=0.15)


def test_he_so_chuyen_loai1_sang_loai2(lasso_res) -> None:
    """eps2 = eps1 * sqrt((lambda*rho + 1)/(lambda*rho)) tai moi buoc."""
    expected = np.sqrt((lasso.LAMBDA * lasso.RHO + 1.0) / (lasso.LAMBDA * lasso.RHO))
    for name in ("S1", "S3", "S5"):
        eps1 = lasso_res["R"][name]["eps1"]
        eps2 = lasso_res["R"][name]["eps2"]
        assert np.allclose(eps2, eps1 * expected, rtol=1e-12)
    assert expected > 1.0  # chung chi loai 2 luon "dat" hon loai 1


def test_lich_trinh_cang_nhieu_lap_noi_thi_sai_so_cang_nho(lasso_res) -> None:
    """T_k=1 (cuc tho) phai cho eps_1 lon hon han T_k=T_ref (gan chinh xac)."""
    eps_tho = lasso_res["R"]["S1"]["eps1"]
    eps_min = lasso_res["R"]["S5"]["eps1"]
    assert np.median(eps_tho) > 100 * np.median(eps_min)


def test_lich_trinh_tot_ve_dich_sat_hon_lich_trinh_tho(lasso_res) -> None:
    """Tang T_k theo k phai thu hep khoang cach toi F* it nhat vai bac."""
    f_star = lasso_res["f_star"]
    gap = {n: lasso_res["R"][n]["Fvals"][-1] - f_star for n in ("S1", "S2", "S3", "S5")}
    assert gap["S1"] > 1e-3          # T_k=1 gan nhu khong tien trien
    assert gap["S2"] < gap["S1"]
    assert gap["S3"] < 1e-6          # T_k ~ log k da du tot
    assert gap["S5"] < 1e-6


def test_gia_tri_toi_uu_tham_chieu_thap_hon_moi_quy_dao(lasso_res) -> None:
    """F* la moc tham chieu, khong quy dao nao duoc xuong duoi qua nhieu."""
    f_star = lasso_res["f_star"]
    for name in lasso_res["R"]:
        assert np.min(lasso_res["R"][name]["Fvals"]) >= f_star - 1e-9


# ===========================================================================
# 4. Khong troi so: ket qua da commit phai tai lap duoc
# ===========================================================================


def test_khop_bang_ket_qua_da_commit(quartic_res) -> None:
    """results/quartic-summary.csv phai khop voi lan chay moi.

    Bat duoc truong hop ai do sua code ma quen chay lai `make figures`,
    khien bang so trong essay/sections/s7 khong con dung voi ma nguon.
    """
    import pandas as pd

    path = quartic.OUT / "quartic-summary.csv"
    if not path.exists():
        pytest.skip("chua co results/quartic-summary.csv")
    df = pd.read_csv(path).set_index("method")

    assert df.loc["IAPPA1 q=7/4", "delta_log_slope_500_2000"] == pytest.approx(
        quartic_res["slope1"], rel=1e-6
    )
    assert df.loc["IAPPA2 q=7/4", "delta_log_slope_500_2000"] == pytest.approx(
        quartic_res["slope2"], rel=1e-6
    )
    assert df.loc["IAPPA1 q=7/4", "final_objective"] == pytest.approx(
        quartic_res["iappa1"][-1], rel=1e-6
    )
    assert df.loc["IAPPA2 q=7/4", "final_objective"] == pytest.approx(
        quartic_res["iappa2"][-1], rel=1e-6
    )
