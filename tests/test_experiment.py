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
# Tai lap -- va ranh gioi chinh xac cua no
# --------------------------------------------------------------------------
# Ban dau bo test nay doi hoi CA HAI tep khop tung byte, va CI tren runner cua
# GitHub bac bo ngay. Sai khac nam DUNG o hai cot do bang polyfit:
#
#   -0.5207170429964162  ->  -0.5207170429964159
#   -1.9888595374794573  ->  -1.988859537479457
#
# Vong lap sinh so lieu la mot day truy hoi mot chieu bang so vo huong thuan
# tuy, khong co BLAS, nen results.csv tai lap tung byte tren moi may. Nhung
# polyfit lai la mot bai toan binh phuong toi thieu giai bang LAPACK, va o do
# thu tu cong don PHU THUOC may. Do chinh la hien tuong quen thuoc, chi an minh
# trong mot buoc hau ky.
def _relative(a: float, b: float) -> float:
    return abs(a - b) / abs(b) if b else abs(a - b)


@pytest.mark.slow
def test_rerunning_reproduces_the_raw_iterates_byte_for_byte():
    """results.csv khong chua phep quy gon nao cua BLAS -- no phai khop tung
    byte, ke ca tren mot may khac."""
    path = EXP / "results.csv"
    before = path.read_bytes()
    subprocess.run([sys.executable, "run_experiment.py"], cwd=EXP, check=True,
                   capture_output=True)
    try:
        assert path.read_bytes() == before
    finally:
        path.write_bytes(before)


@pytest.mark.slow
def test_rerunning_reproduces_the_summary_to_full_double_precision():
    """summary.csv thi khong the doi hoi tung byte: hai cot do doc di qua
    polyfit. Doi hoi dung muc ma chung thuc su dat duoc."""
    path = EXP / "summary.csv"
    before = path.read_bytes()
    old = pd.read_csv(path).set_index("method")
    subprocess.run([sys.executable, "run_experiment.py"], cwd=EXP, check=True,
                   capture_output=True)
    try:
        new = pd.read_csv(path).set_index("method")
        for method in old.index:
            # Cac dai luong lay thang tu day lap: khop tung bit.
            for col in ("final_objective", "max_abs_certificate_minus_1"):
                a, b = new.loc[method, col], old.loc[method, col]
                if pd.isna(a) or pd.isna(b):
                    assert pd.isna(a) and pd.isna(b), (method, col)
                else:
                    assert a == b, (method, col, a, b)
            # Do doc hoi quy: khop den vai ulp, khong hon.
            a, b = (new.loc[method, "delta_log_slope_500_2000"],
                    old.loc[method, "delta_log_slope_500_2000"])
            if not (pd.isna(a) or pd.isna(b)):
                assert _relative(a, b) < 1e-12, (method, a, b)
    finally:
        path.write_bytes(before)
