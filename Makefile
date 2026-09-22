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
#  make test       kiểm chứng số học (các khẳng định toán học phải còn đúng)
#  make lint       chạy ruff trên mã nguồn Python
#  make verify     kiểm chứng chặn hội tụ của báo cáo bằng Octave (không cần MATLAB)
#  make matexp     chạy lại thực nghiệm MATLAB/Octave (Octave: bỏ qua phần vẽ hình)
#  make clean      xoá file trung gian, giữ lại PDF
#  make distclean  xoá cả PDF
# ============================================================================

LATEXMK   ?= latexmk
LATEXMKFLAGS ?= -pdf -interaction=nonstopmode -halt-on-error
PYTHON    ?= python3
OCTAVE    ?= octave --no-gui --quiet

.PHONY: all report slides transcript essay notes figures matexp test verify lint format clean distclean help

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
figures:
	$(PYTHON) code/python/quartic_iappa.py
	$(PYTHON) code/python/lasso_iappa.py

# Thuc nghiem MATLAB. Phan TINH TOAN chay duoc bang GNU Octave (khong can
# giay phep MATLAB) -- do la phan CI kiem. Phan VE HINH thi Octave khong chay
# duoc tren cac may da thu, nen script tu bo qua; ba tep figures/*.pdf cua bao
# cao van phai sinh bang MATLAB.
matexp:
	cd code/matlab && $(OCTAVE) --eval "main_experiment"

# Kiem chung so hoc: khac `make figures` o cho no CO THE THAT BAI.
# Kiem cac dang thuc dong cua delta_k, hai chan hoi tu, chung chi sai so,
# va tinh dung dan cua prox / soft-threshold / lien hop Fenchel.
test:
	cd code/python && $(PYTHON) -m pytest -q

# Khang dinh trung tam cua bao cao la mot chan hoi tu. Muc tieu nay chay thuat
# toan va kiem tra chan do KHONG bi vi pham, tai moi buoc lap ngoai, voi moi
# lich trinh lap noi -- va kiem ca gia thiet cua chinh dinh ly. Chay khong can
# giao dien do hoa, nen dung duoc truc tiep trong CI.
verify:
	cd code/matlab && $(OCTAVE) --eval "exit(~verify_bounds())"

lint:
	$(PYTHON) -m ruff check .

format:
	$(PYTHON) -m ruff check --fix .

clean:
	$(LATEXMK) -c main.tex slide.tex transcript.tex
	cd essay && $(LATEXMK) -c essay.tex
	cd docs  && $(LATEXMK) -c ghi-chu-bai-bao.tex

distclean:
	$(LATEXMK) -C main.tex slide.tex transcript.tex
	cd essay && $(LATEXMK) -C essay.tex
	cd docs  && $(LATEXMK) -C ghi-chu-bai-bao.tex

help:
	@sed -n '2,17p' $(MAKEFILE_LIST)
