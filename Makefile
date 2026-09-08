# ============================================================================
#  Inexact and Accelerated Proximal Point Algorithms — Nhóm 10
#
#  make            biên dịch toàn bộ 4 tài liệu PDF
#  make report     chỉ báo cáo chính (main.pdf)
#  make slides     chỉ slide (slide.pdf)
#  make transcript chỉ kịch bản thuyết trình (transcript.pdf)
#  make essay      chỉ tiểu luận rút gọn (essay/essay.pdf)
#  make notes      ghi chú đọc bài báo (docs/ghi-chu-bai-bao.pdf)
#  make figures    chạy lại thực nghiệm Python, sinh lại hình + bảng số liệu
#  make clean      xoá file trung gian, giữ lại PDF
#  make distclean  xoá cả PDF
# ============================================================================

LATEXMK   ?= latexmk
LATEXMKFLAGS ?= -pdf -interaction=nonstopmode -halt-on-error
PYTHON    ?= python3

.PHONY: all report slides transcript essay notes figures clean distclean help

all: report slides transcript essay

report:     main.pdf
slides:     slide.pdf
transcript: transcript.pdf

main.pdf slide.pdf transcript.pdf: %.pdf: %.tex metadata.tex refs.bib $(wildcard content/*.tex)
	$(LATEXMK) $(LATEXMKFLAGS) $<

essay: essay/essay.pdf
essay/essay.pdf: essay/essay.tex essay/preamble.tex essay/refs.bib $(wildcard essay/sections/*.tex)
	cd essay && $(LATEXMK) $(LATEXMKFLAGS) essay.tex

notes: docs/ghi-chu-bai-bao.pdf
docs/ghi-chu-bai-bao.pdf: docs/ghi-chu-bai-bao.tex
	cd docs && $(LATEXMK) $(LATEXMKFLAGS) ghi-chu-bai-bao.tex

# Sinh lại hình + bảng số liệu từ thực nghiệm Python.
# Hình MATLAB (figures/convergence-comparison.pdf, ...) phải chạy riêng bằng
# MATLAB: cd code/matlab && main_experiment
figures:
	$(PYTHON) code/python/quartic_iappa.py
	$(PYTHON) code/python/lasso_iappa.py

clean:
	$(LATEXMK) -c main.tex slide.tex transcript.tex
	cd essay && $(LATEXMK) -c essay.tex
	cd docs  && $(LATEXMK) -c ghi-chu-bai-bao.tex

distclean:
	$(LATEXMK) -C main.tex slide.tex transcript.tex
	cd essay && $(LATEXMK) -C essay.tex
	cd docs  && $(LATEXMK) -C ghi-chu-bai-bao.tex

help:
	@sed -n '2,14p' $(MAKEFILE_LIST)
