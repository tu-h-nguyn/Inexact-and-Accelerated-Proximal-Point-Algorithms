"""Kiem chung cac khang dinh ly thuyet tu ket qua da cong bo trong test_01/.

Khac voi phan MATLAB (xem matlab/verify_bounds.m), thi nghiem Python nay tai
lap duoc TUNG BIT: no chi dung numpy tren mot bai toan mot chieu tat dinh,
khong co bo sinh so ngau nhien va khong co phep thu gon song song nao lam doi
thu tu cong don. Vi vay o day co the vua kiem tra khang dinh, vua kiem tra
chinh su tai lap -- dieu ma phia MATLAB khong the doi hoi.
"""
from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / "test_01"


@pytest.fixture(scope="module")
def summary() -> pd.DataFrame:
    return pd.read_csv(EXP / "summary.csv").set_index("method")


@pytest.fixture(scope="module")
def results() -> pd.DataFrame:
    return pd.read_csv(EXP / "results.csv")


# --------------------------------------------------------------------------
# Bac hoi tu: khang dinh trung tam cua bao cao
# --------------------------------------------------------------------------
def test_type2_error_decays_at_the_predicted_order(summary):
    """Sai so loai 2 phai suy giam bac 2 -- do la ly do no "khoi phuc" duoc
    toc do cua thuat toan chinh xac."""
    slope = summary.loc["IAPPA2 q=7/4", "delta_log_slope_500_2000"]
    assert slope == pytest.approx(-2.0, abs=0.05), slope


def test_type1_error_decays_much_more_slowly(summary):
    """Voi cung mot lich trinh sai so q = 7/4, loai 1 chi dat khoang k^{-1/2}.
    Day chinh la khoang cach ma Chuong 4 phan tich."""
    slope = summary.loc["IAPPA1 q=7/4", "delta_log_slope_500_2000"]
    assert slope == pytest.approx(-0.5, abs=0.05), slope


def test_type2_is_at_least_three_orders_faster(summary):
    s1 = summary.loc["IAPPA1 q=7/4", "delta_log_slope_500_2000"]
    s2 = summary.loc["IAPPA2 q=7/4", "delta_log_slope_500_2000"]
    assert s2 < s1 - 1.0, (s1, s2)


# --------------------------------------------------------------------------
# Chung chi: dai luong phai bang 1 neu cac dang thuc duoc dung
# --------------------------------------------------------------------------
@pytest.mark.parametrize("method", ["IAPPA1 q=7/4", "IAPPA2 q=7/4"])
def test_certificate_is_one(summary, method):
    """`cert` la mot dang thuc phai thoa dung bang 1 tai moi vong lap neu cac
    buoc bien doi trong chung minh la dung. Lech khoi 1 la loi cai dat."""
    dev = summary.loc[method, "max_abs_certificate_minus_1"]
    assert dev < 1e-8, dev


# --------------------------------------------------------------------------
# Tinh nhat quan cua ban than day lap
# --------------------------------------------------------------------------
def test_all_objective_gaps_are_nonnegative(results):
    """F(x_k) - F_* khong the am: F_* la gia tri toi tieu."""
    for col in ("exact", "IAPPA1_q1.75", "IAPPA2_q1.75"):
        assert (results[col] >= 0).all(), col


def test_error_components_are_nonnegative_and_finite(results):
    for col in ("delta1_bound_component", "delta2_bound_component"):
        v = results[col]
        assert v.notna().all() and (v >= 0).all() and v.map(math.isfinite).all(), col


def test_exact_algorithm_is_the_most_accurate_of_the_three(results):
    """Thuat toan chinh xac khong the thua ca hai ban khong chinh xac: neu no
    thua, hoac cai dat sai, hoac "chinh xac" khong thuc su chinh xac."""
    tail = results.tail(200)
    assert tail["exact"].median() <= tail["IAPPA1_q1.75"].median()


def test_final_objective_ordering_matches_the_report(summary):
    exact = summary.loc["Exact accelerated", "final_objective"]
    t1 = summary.loc["IAPPA1 q=7/4", "final_objective"]
    t2 = summary.loc["IAPPA2 q=7/4", "final_objective"]
    # Loai 1 la ban te nhat sau 2000 vong; loai 2 bat kip thuat toan chinh xac.
    assert t1 > exact and t1 > t2, (exact, t1, t2)


# --------------------------------------------------------------------------
# Tai lap tung bit -- chay lai that su
# --------------------------------------------------------------------------
@pytest.mark.slow
def test_rerunning_the_experiment_reproduces_the_committed_csvs():
    before = {p: p.read_bytes() for p in (EXP / "results.csv", EXP / "summary.csv")}
    subprocess.run([sys.executable, "run_experiment.py"], cwd=EXP, check=True,
                   capture_output=True)
    try:
        for p, old in before.items():
            assert p.read_bytes() == old, f"{p.name} thay doi sau khi chay lai"
    finally:
        for p, old in before.items():
            p.write_bytes(old)
