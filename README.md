# Thuật toán điểm gần kề tăng tốc có sai số

**Inexact and Accelerated Proximal Point Algorithms** — báo cáo, slide, kịch bản thuyết trình và thực nghiệm số tái lập được, trình bày lại công trình của **Saverio Salzo & Silvia Villa** (*Journal of Convex Analysis* 19, 2012).

[![build](https://github.com/tu-h-nguyn/inexact-and-accelerated-proximal-point-algorithms/actions/workflows/build.yml/badge.svg)](https://github.com/tu-h-nguyn/inexact-and-accelerated-proximal-point-algorithms/actions/workflows/build.yml)
![LaTeX](https://img.shields.io/badge/LaTeX-pdflatex-008080)
![Python](https://img.shields.io/badge/Python-3.11-3776ab)
![MATLAB](https://img.shields.io/badge/MATLAB-R2020b%2B-e16737)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **In English.** This repository reproduces and extends the analysis of Salzo & Villa (2012) on
> inexact accelerated proximal point algorithms. The paper identifies a subtle flaw in Güler's 1992
> convergence proof for the inexact case and rebuilds the analysis on a flexible *estimate sequence*
> framework. The deliverables are a 34-page report, a 45-slide deck, a presentation script, a
> condensed standalone essay, and two fully reproducible numerical experiments (Python + MATLAB)
> that verify both theoretical bounds numerically. All documents and figures are rebuilt in CI.

---

## Câu hỏi trung tâm

Toán tử gần kề

$$\operatorname{prox}_{\lambda F}(y) = \arg\min_x \Big\{ F(x) + \tfrac{1}{2\lambda}\lVert x-y\rVert^2 \Big\}$$

hiếm khi có công thức đóng: trong thực tế mỗi bước lặp phải **giải một bài toán con bằng một thuật toán khác**, tức là chỉ thu được một xấp xỉ. Năm 1992 Güler kết hợp PPA với ngoại suy Nesterov và công bố tốc độ $\mathcal{O}(1/k^2)$ — kể cả khi tính xấp xỉ.

Salzo & Villa chỉ ra rằng chứng minh đó có một lỗ hổng tinh tế: lập luận **ngầm sử dụng dưới gradient tại điểm gần kề chính xác** — một đại lượng không thể tính được khi ta chỉ có xấp xỉ. Bài báo xây dựng lại toàn bộ phân tích trên khung *dãy ước lượng* và trả lời: **tốc độ hội tụ còn giữ được hay không phụ thuộc vào việc bộ giải bài toán con cấp cho ta loại chứng nhận sai số nào.**

Với $\varepsilon_k=\mathcal{O}(1/k^q)$, hai thuật toán cho hai bức tranh hoàn toàn khác nhau:

| Loại sai số | Điều kiện | Tốc độ tốt nhất đạt được | Ý nghĩa thực tế |
|---|---|---|---|
| **Loại 1** (AT1) → IAPPA1 | $0 \in \partial_{\varepsilon^2/2\lambda}\Phi_\lambda(z)$ | $\mathcal{O}(1/k)$ — và chỉ khi $q>2$; $q<2$ cho $\mathcal{O}(1/k^{2q-3})$ | Tiêu chí dừng tự nhiên nhất, nhưng **mất sạch lợi thế tăng tốc**: dù $\varepsilon_k$ giảm nhanh đến đâu cũng không vượt được $\mathcal{O}(1/k)$ |
| **Loại 2** (AT2) → IAPPA2 | $\frac{y-z}{\lambda} \in \partial_{\varepsilon^2/2\lambda}F(z)$ | $\mathcal{O}(1/k^2)$ khi $q>3/2$; $(x_k)$ vẫn tối tiểu hóa với mọi $q>1/2$ | Khắt khe hơn, nhưng **khôi phục đầy đủ tốc độ bậc hai** |
| **Loại 3** (AT3) | $d(0,\partial\Phi_\lambda(z)) \le \varepsilon/\lambda$ | — | Tương đương prox chính xác của một đầu vào bị nhiễu |

Kết luận đắt giá nhất: **tăng tốc không miễn phí.** Nếu chỉ đo được sai số theo loại 1, thuật toán tăng tốc chạy không nhanh hơn PPA thường; muốn giữ $\mathcal{O}(1/k^2)$ thì bộ giải bài toán con phải cấp được chứng nhận loại 2.

---

## Kết quả thực nghiệm

### 1. Kiểm chứng bậc hội tụ trên hàm bậc bốn (Python, $\mathcal{H}=\mathbb{R}$)

$F(x)=x^4/4$, sai số dựng **đúng trên biên** ngân sách $\varepsilon_k = 0{,}2/(k+1)^{7/4}$ để không vô tình tạo nhiễu có lợi. Hồi quy $\log\delta_k$ theo $\log k$ trên đoạn $500 \le k \le 2000$:

| Phương pháp | $F(x_N)-F_\star$ | Hệ số góc của $\delta_k$ | Bậc lý thuyết | Sai lệch chứng chỉ $\lvert\rho_k-1\rvert$ |
|---|---|---|---|---|
| Tham chiếu (giải gần chính xác) | $5{,}266\times10^{-13}$ | — | — | — |
| **IAPPA1**, $q=7/4$ | $6{,}067\times10^{-10}$ | $-0{,}521$ | $-1/2$ | $8{,}0\times10^{-12}$ |
| **IAPPA2**, $q=7/4$ | $7{,}440\times10^{-16}$ | $-1{,}989$ | $-2$ | $7{,}8\times10^{-10}$ |

Hai hệ số góc đo được khớp với hai bậc lý thuyết đến chữ số thứ hai — chênh lệch giữa loại 1 và loại 2 hiện ra rõ ràng bằng số.

<p align="center">
  <img src="figures/quartic-error-accumulation.png" width="70%" alt="Tích lũy sai số của IAPPA1 và IAPPA2">
</p>

### 2. Bài toán elastic-net / LASSO không có prox đóng (MATLAB + Python)

$F(x)=\tfrac12\lVert Ax-b\rVert^2 + \tfrac{\rho}{2}\lVert x\rVert^2 + \mu\lVert x\rVert_1$ với $n=80$, $m=30$, độ thưa $s=8$. Ở đây $\operatorname{prox}_{\lambda F}$ **thật sự không có công thức đóng**, nên bài toán con được giải bằng FISTA nội — số vòng lặp nội $T_k$ chính là "nút vặn" sinh sai số có kiểm soát.

Điểm đáng chú ý về mặt cài đặt: **chỉ cần một lần giải gần đúng** là thu được đồng thời cả hai chứng nhận, vì $F$ lồi mạnh hệ số $\rho$ cho phép chuyển thẳng

$$\varepsilon^{(1)} = \sqrt{2\lambda\,\delta_{\text{inner}}}, \qquad \varepsilon^{(2)} = \varepsilon^{(1)}\sqrt{\tfrac{\lambda\rho+1}{\lambda\rho}}.$$

Nhờ vậy cùng một quỹ đạo $(x_k)$ đi qua được cả hai bộ sổ sách kế toán sai số và hai chặn lý thuyết được so sánh trên **cùng một dữ liệu**.

<p align="center">
  <img src="figures/lasso-convergence-comparison.png" width="49%" alt="So sánh tốc độ hội tụ theo lịch trình lặp nội">
  <img src="figures/lasso-bound-verification.png" width="49%" alt="Kiểm chứng số học hai chặn lý thuyết">
</p>

Kết quả: lịch trình $T_k \sim \sqrt{k}$ giữ được tốc độ $\mathcal{O}(1/k^2)$ gần như ngang với việc giải chính xác, trong khi $T_k$ cố định ở mức thấp tụt về $\mathcal{O}(1/k)$. Cả hai chặn lý thuyết đều được kiểm chứng: **0/121 bước vi phạm** ở cả loại 1 lẫn loại 2.

---

## Cấu trúc repository

```
├── main.tex              Báo cáo đầy đủ (34 trang) — điều phối chương/mục, \input từ content/
├── slide.tex             Slide Beamer (45 trang) — dùng lại đúng các file trong content/
├── transcript.tex        Kịch bản thuyết trình dạng screenplay (21 trang)
├── metadata.tex          Tên đề tài, nhóm, thành viên, GVHD, ngày — sửa DUY NHẤT ở đây
├── refs.bib              Thư mục tài liệu tham khảo (BibTeX)
├── content/              Thân bài dùng chung cho cả ba tài liệu trên (mọi file ở đây đều được biên dịch)
│
├── essay/                Tiểu luận rút gọn, độc lập (12 trang, 8 mục, bib riêng)
├── docs/                 Tài liệu phụ trợ, KHÔNG thuộc bản nộp (xem docs/README.md)
│
├── code/
│   ├── matlab/           Thực nghiệm LASSO (MATLAB hoặc Octave) + verify_bounds.m (38 phép kiểm)
│   └── python/           quartic_iappa.py + lasso_iappa.py + test_iappa.py (52 phép kiểm)
│
├── figures/              Toàn bộ hình dùng trong tài liệu (.pdf vector + .png)
├── results/              Bảng số liệu .csv do thực nghiệm sinh ra
│
├── Makefile              make all / test / verify / figures / matfigures / clean
└── .github/workflows/    CI: build 5 PDF + kiểm chứng số học (Python + Octave) + chạy lại thực nghiệm
```

---

## Biên dịch

```bash
make all        # main.pdf, slide.pdf, transcript.pdf, essay/essay.pdf
make notes      # docs/ghi-chu-bai-bao.pdf
make report     # chỉ báo cáo
make clean      # xoá file trung gian, giữ PDF
```

Không có `make`? Chạy tay (lặp `pdflatex` để mục lục và `\cite` ổn định):

```bash
pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
```

**Yêu cầu:** TeX Live với `texlive-lang-other` (gói `vietnam`/`vntex` cho tiếng Việt), `texlive-latex-extra`, `texlive-science`, `latexmk`. Trên Overleaf: upload cả repo rồi đặt `main.tex` làm main document.

## Tái lập thực nghiệm số

```bash
pip install -r code/python/requirements.txt
make test           # ~25 giây, 52 phép kiểm chứng số học (Python)
make figures        # ~13 giây, ghi đè figures/*.png và results/*.csv
make verify         # ~40 giây, 38 phép kiểm chặn hội tụ (Octave, không cần MATLAB)
make matexp         # chạy lại thực nghiệm MATLAB (Octave chạy được phần tính toán)
```

`make test` **không** kiểm tra "script có chạy không" — việc đó `make figures` đã làm.
Nó kiểm tra **các khẳng định toán học còn đúng hay không**, và nó *thất bại được*:

| Nhóm | Kiểm cái gì |
|---|---|
| Khối cơ bản | `prox`, `soft_threshold` thoả điều kiện tối ưu bậc một; `F*` đúng là liên hợp Fenchel (đẳng thức Fenchel–Young đạt tại $v=F'(x)$) |
| Chứng chỉ sai số | Điểm dựng nằm **đúng** trên biên ngân sách ($\rho_k=1$), không tự cho mình xấp xỉ tốt hơn mức khai báo |
| Dạng đóng của $\delta_k$ | $\delta_k=\frac{\beta_k}{2}\sum\frac{\varepsilon_i^2}{\lambda_i\beta_{i+1}}$ (IAPPA2) và $\delta_k=\frac{A\beta_k}{2}\sum\eta_{i+1}^2$ (IAPPA1) |
| Chặn hội tụ | Không bước nào vi phạm chặn loại 1 hoặc loại 2 |
| Bậc hội tụ | Hệ số góc của $\delta_k$ khớp $-1/2$ và $-2$ trong sai số $0{,}1$ |
| Không trôi số | `results/quartic-summary.csv` đã commit khớp với lần chạy mới |

Bộ kiểm này đã được **kiểm tra ngược bằng đột biến**: cố tình gài 10 lỗi vào mã nguồn
(đổi hệ số trong công thức cập nhật $\delta$, bỏ nhân $(1-\alpha)$, sai hệ số chuyển
loại 1 → loại 2, dùng nửa ngân sách sai số để kết quả đẹp giả, chọn điểm biên có lợi
thay vì bất lợi…) — **cả 10 đều bị bắt**.

### Bản MATLAB — và vì sao nó được kiểm khác đi

Bản MATLAB sinh 3 hình vector `.pdf` dùng trong báo cáo. **Phần tính toán** chạy
được bằng GNU Octave, không cần giấy phép MATLAB — và CI kiểm đúng điều đó:

```bash
cd code/matlab
octave --no-gui --quiet --eval "main_experiment"     # hoặc: make matexp
```

**Phần vẽ hình thì không.** `ft_text_renderer` của Octave không nạp được font trên
cả hai môi trường đã thử (kể cả sau khi cài `gnuplot-nox` và chỉ định font DejaVu),
nên script tự phát hiện và bỏ qua phần đồ thị thay vì chết giữa chừng — kết quả số
vẫn in ra đầy đủ. Ba tệp `figures/*.pdf` của báo cáo vẫn phải sinh bằng MATLAB.
`export_fig_pdf.m` chọn `exportgraphics` hay `print` tuỳ môi trường, nên đường vẽ
hình sẽ tự chạy nếu Octave trên máy bạn không dính lỗi font này.

Nhưng nó **không thể** bị kiểm bằng cách so từng chữ số như bản Python. MATLAB và
Octave dùng hai bộ sinh số ngẫu nhiên khác nhau, nên cùng một hạt giống vẫn cho hai
*thực thể bài toán* khác nhau; so con số sẽ đỏ ngay khi đổi môi trường mà không bắt
được lỗi thật nào. Vì vậy `make verify` kiểm những điều phải đúng trên **mọi** thực
thể — tức chính là nội dung của các định lý:

| Nhóm | Kiểm cái gì |
|---|---|
| Chặn hội tụ (Định lý 3.2) | **0 vi phạm / 121 bước ngoài**, cho cả 5 lịch trình lặp nội |
| Giả thiết của chính định lý | Định lý đòi $a \le \alpha_k^2/\bigl((1-\alpha_k)A_k\lambda_k\bigr) \le 2$; công thức đóng trong `run_outer.m` làm tỉ số này **đúng bằng 1** — kiểm tới $7\times10^{-16}$ |
| Đẳng thức đóng | $\varepsilon_2/\varepsilon_1=\sqrt{(\lambda\rho+1)/(\lambda\rho)}$; truy hồi $A_{k+1}=(1-\alpha_k)A_k$; tỉ số $\delta_2/\delta_1$ tại bước 1 **đúng bằng 26** |
| Điểm giao loại 1 / loại 2 | Loại 2 **thiệt thòi 26 lần** ở bước đầu, chỉ vượt lên từ bước 9 — "khôi phục bậc hội tụ" là phát biểu *tiệm cận*, và bộ kiểm khẳng định đúng điều đó |
| Tác dụng của độ chính xác nội | $T_k=1$ để lại khoảng cách $1{,}0\times10^{-1}$, gần bảy bậc tệ hơn $1{,}5\times10^{-8}$; lịch trình *tăng dần* lấy lại gần hết |

Bộ kiểm này **độc lập với hạt giống** (đạt 38/38 với hạt 1, 7 và 42) và đã được kiểm
ngược bằng đột biến: đổi hệ số 4 thành 3,5 trong công thức $\alpha_k$ thì bị bắt ngay.

Hai bản cho cùng kết luận định tính; con số tuyệt đối lệch nhau vì lý do trên. Cả hai
đều chạy trong CI, nên **người đọc không có MATLAB vẫn tái lập được toàn bộ kết quả**.

---

## Quy ước biên soạn

### Một nội dung, ba đầu ra

`main.tex`, `slide.tex`, `transcript.tex` cùng `\input` các file trong `content/`; mỗi bản chỉ hiện phần dành cho mình nhờ package `comment`:

| Environment | Báo cáo | Slide | Kịch bản |
|---|:---:|:---:|:---:|
| `essayonly` | ✅ | ❌ | ❌ |
| `slidesonly` | ❌ | ✅ | ❌ |
| `scriptonly` | ❌ | ❌ | ✅ |

Chứng minh chi tiết và lập luận dài nằm trong `essayonly` (chỉ vào báo cáo); slide chỉ giữ phát biểu định lý và công thức chính. **Muốn biết chính xác nội dung nào đang hiện trên slide** thì đọc phần *không* nằm trong `\begin{essayonly}...\end{essayonly}`.

### Thêm một mục mới

1. Tạo file trong `content/` theo quy ước `chap{chương}_{thứ tự}_{tên gợi nhớ}.tex`, chỉ viết thân bài (**không** viết `\chapter`/`\section` trong đó).
2. Thêm `\section{...}` rồi `\input{content/ten_file.tex}` — làm ở **cả `main.tex` lẫn `slide.tex`**.
3. Cần `\ref{}` từ nơi khác thì gắn `\label{...}` ngay sau `\section`.

### Quy tắc chung

- Overleaf chỉ cho 2 người sửa cùng lúc nên nhóm dùng repo này thay thế — **commit + push thường xuyên**.
- Không commit file build (`.aux`, `.log`, `.bbl`, PDF sinh ra…) — đã có trong `.gitignore`. Hình trong `figures/` là dữ liệu đầu vào nên vẫn được commit.
- Đánh nhãn công thức bằng tên gợi nhớ (`\autotag{eq:ten-nhan}`), **không** dựa vào số thứ tự — bộ đếm equation khác nhau giữa báo cáo và slide.

---

## Trạng thái

| Tài liệu | Số trang | Trạng thái |
|---|---:|---|
| `main.pdf` — Báo cáo | 34 | ✅ Sạch: 0 tham chiếu treo, 0 chỗ tràn lề |
| `slide.pdf` — Slide | 45 | ✅ Biên dịch sạch |
| `transcript.pdf` — Kịch bản | 21 | ✅ Biên dịch sạch, canh lề trái theo quy ước kịch bản |
| `essay/essay.pdf` — Tiểu luận | 12 | ✅ Biên dịch sạch, số liệu khớp `results/` |
| `docs/ghi-chu-bai-bao.pdf` | 9 | ✅ Biên dịch sạch |

Cả năm tài liệu: **0 tham chiếu/trích dẫn treo, 0 chỗ tràn lề quá 20pt**. CI kiểm tra lại toàn bộ điều trên ở mỗi lần push, chạy lại hai thực nghiệm số và **chạy 52 phép kiểm chứng toán học** — nếu một kết luận trong báo cáo không còn đúng với mã nguồn, CI đỏ.

---

## Nhóm thực hiện

**Nhóm 10** — Học phần Phương pháp số trong Tối ưu, Khoa Toán – Tin học, Trường Đại học Khoa học Tự nhiên, ĐHQG-HCM.

| Thành viên | MSSV |
|---|---|
| Nguyễn Hoàng Tú | 23110220 |
| Bùi Công Hoàng Vũ | 23110223 |
| Huỳnh Trung Kiên | 22110091 |
| Nguyễn Ngọc Diễm Quỳnh | 21110167 |

Giảng viên hướng dẫn: **TS. Nguyễn Đăng Khoa**

## Tài liệu gốc

> Saverio Salzo, Silvia Villa. *Inexact and Accelerated Proximal Point Algorithms.*
> Journal of Convex Analysis **19** (2012), no. 4, 1167–1192.

Các tài liệu nền tảng khác (Güler 1992, Nesterov 2004, Beck–Teboulle 2009, Rockafellar 1976…) có đầy đủ trong [`refs.bib`](refs.bib).

## Giấy phép

Mã nguồn và phần văn bản do nhóm biên soạn: [MIT](LICENSE). Nội dung toán học được trình bày lại từ bài báo gốc — bản quyền thuộc về các tác giả và nhà xuất bản.
