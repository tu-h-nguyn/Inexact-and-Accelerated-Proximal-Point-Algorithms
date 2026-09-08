# Report-Slide-Numericals-for-Optimizations-Topic-10

Repo chứa báo cáo (`main.tex`), slide (`slide.tex`) và kịch bản thuyết trình (`transcript.tex`) cho đề tài, cùng mã nguồn MATLAB minh họa số (`matlab/`). Ghi chú dưới đây để mọi người nắm được luồng làm việc, không cần hỏi lại.

## Quy tắc làm việc chung
- Overleaf chỉ cho tối đa 2 người edit cùng lúc, nên cả nhóm dùng repo này thay thế. Vì vậy **nhớ commit + push thường xuyên** để tránh conflict và mất bài.
- Repo dùng GitHub nên mỗi người tự chuẩn bị compiler riêng (TeXworks, Overleaf tải project lên, VS Code + LaTeX Workshop, v.v.) — không có server build chung.
- Không commit file build (`.pdf`, `.aux`, `.log`, `.bbl`,...) — các file này đã có trong `.gitignore`, chỉ commit source (`.tex`, `.bib`, hình ảnh). Ngoại lệ: 3 file `.pdf` trong `matlab/` là hình xuất từ script, được `\includegraphics` thẳng vào báo cáo nên có commit.

## Cấu trúc thư mục
```
main.tex           <- báo cáo: khai báo chương/mục, rồi \input nội dung từ content/
slide.tex          <- slide Beamer, dùng lại đúng các file trong content/ (không tự chứa nội dung riêng)
transcript.tex      <- kịch bản thuyết trình dạng screenplay (scene/action/character/speech), đọc song song với slide
metadata.tex       <- tên đề tài, tên nhóm, thành viên, GVHD, ngày tháng (sửa ở đây, không sửa trong main.tex/slide.tex)
refs.bib           <- danh sách tài liệu tham khảo (BibTeX), cite bằng \citep{key} / \cite{key}
content/           <- mỗi file là nội dung của MỘT mục hoặc tiểu mục, không tự chứa \chapter/\section
matlab/            <- code MATLAB thực nghiệm số (xem mục riêng bên dưới)
```

## Luồng biên soạn: main.tex điều phối, content/ chứa nội dung
Tiêu đề chương (`\chapter`), mục (`\section`), tiểu mục (`\subsection`) được khai báo **trong `main.tex`** (và tương ứng trong `slide.tex`). Mỗi file trong `content/` chỉ chứa phần thân bài (định lý, chứng minh, hình vẽ...) của đúng một mục/tiểu mục đó, và được chèn vào bằng `\input{content/ten_file.tex}`.

Quy ước đặt tên file: `chap{số chương}_{số thứ tự}_{tên gợi nhớ}.tex`, ví dụ `chap3_2_xay_dung.tex` là file thứ 2 của Chương 3.

**Muốn thêm một mục/tiểu mục mới:**
1. Tạo file mới trong `content/` theo quy ước tên trên, chỉ viết nội dung (không viết `\chapter`/`\section` trong đó).
2. Thêm dòng `\section{Tên mục}` (hoặc `\subsection{...}`) ngay trước, rồi `\input{content/ten_file_vua_tao.tex}` ngay sau — làm ở **cả `main.tex` lẫn `slide.tex`** vì hai file dùng chung nội dung.
3. Nếu mục nào cần được `\ref{}` từ nơi khác, gắn `\label{...}` ngay sau lệnh `\section`/`\subsection` đó.

## Cơ chế lọc nội dung essayonly / slidesonly / scriptonly
Cả ba file (`main.tex`, `slide.tex`, `transcript.tex`) `\input` chung các file trong `content/`, nhưng mỗi phiên bản chỉ hiện một phần nhờ package `comment`:

| Environment      | main.tex (báo cáo) | slide.tex (slide) | transcript.tex (kịch bản) |
|---|---|---|---|
| `essayonly`   | hiện | ẩn | ẩn |
| `slidesonly`  | ẩn | hiện | ẩn |
| `scriptonly`  | ẩn | ẩn | hiện |

Dùng để: chứng minh chi tiết / lập luận dài chỉ để trong `essayonly` (chỉ vào báo cáo), còn slide chỉ giữ phát biểu định lý/công thức chính. Muốn biết **chính xác nội dung nào đang hiện trên slide** thì đọc phần **không** nằm trong `\begin{essayonly}...\end{essayonly}` của file `content/` tương ứng.

## transcript.tex — kịch bản thuyết trình
File screenplay riêng, không `\input` từ `content/` mà viết tay để khớp với đúng nội dung hiện trên slide (xem cơ chế lọc ở trên). Cấu trúc:
- `\scene{}`, `\action{}`, `\character{}`, `\parenthetical{}` — định dạng kịch bản.
- `\begin{speech}...\end{speech}` — từng khối lời thoại, chia nhỏ ~25 từ/khối để dễ đọc khi trình chiếu.
- `\block{N}` — in nhãn "(~N từ)" phía trên mỗi khối, N phải khớp đúng số từ thực tế trong khối (đã được rà soát tự động, không còn khối nào lệch số).
- Mọi khoảng trống nội dung trên slide trước đây (ví dụ: phần Ví dụ phép chiếu, Bổ đề lem:core, Định lý 4.1/4.3 + IAPPA1, Định lý 4.2 + IAPPA2) đã được thuyết minh đầy đủ, bám sát nội dung thật sự hiển thị trên slide.

## matlab/ — thực nghiệm số
`main_experiment.m` là script chính, chạy bài toán LASSO tổng hợp (n=80, m=30, sparsity s=8) qua thuật toán điểm gần kề tăng tốc không chính xác, dùng các hàm phụ trợ:
- `F_obj.m` — hàm mục tiêu (L1 + ridge).
- `soft_threshold.m` — toán tử prox của chuẩn 1.
- `inner_fista_trace.m` — bộ giải FISTA nội, trả về vết lặp.
- `run_outer.m` — vòng lặp ngoài chính (IAPPA), tích lũy `delta1`/`delta2`, `eps1`/`eps2`.
- `run_classical.m` — đường cơ sở PPA cổ điển không tăng tốc, để so sánh.

Script xuất ra 3 hình (đã `\includegraphics` vào `content/outro.tex`, mục "Tổng kết"):
1. **`So sanh toc do hoi tu.pdf`** — so sánh tốc độ hội tụ giữa PPA cổ điển và các lịch trình lặp nội khác nhau (`T_k` hằng số / log / sqrt(k) / chính xác), đối chiếu với `O(1/k)` và `O(1/k^2)`.
2. **`Kiem chung so hoc chan hoi tu cua Dinh ly 3-2.pdf`** — kiểm chứng số học chận hội tụ (Định lý 3.2 / chương 4), so sánh sai số thực tế với hai chận lý thuyết loại 1 và loại 2.
3. **`Do chinh xac loai 1 thuc te dat duoc theo tung lich trinh lap noi.pdf`** — độ suy giảm của `\eps_k` (sai số loại 1) thực tế theo từng lịch trình lặp nội, minh họa điều kiện suy giảm cần thiết để giữ tốc độ hội tụ bậc 2.

Muốn tái tạo hình: chạy `main_experiment.m` trong MATLAB, các file `.pdf` sẽ được ghi đè tại chỗ (cùng thư mục `matlab/`).

## Trạng thái hiện tại
Toàn bộ nội dung đã hoàn chỉnh và khớp nhau giữa 3 file:
- **Báo cáo (`main.tex`)**: đầy đủ — Lời nói đầu, Chương 1 "Cơ sở toán học" (4 mục), Chương 2 "Dãy ước lượng Nesterov" (3 mục), Chương 3 "Bậc hội tụ của thuật toán" (2 mục: sai số loại 1, loại 2), Tổng kết (rút gọn thuật toán, code MATLAB, 3 hình thực nghiệm), tài liệu tham khảo.
- **Slide (`slide.tex`)**: 41 trang, đi theo đúng cấu trúc trên, chỉ giữ phần `slidesonly`.
- **Kịch bản (`transcript.tex`)**: khớp 1-1 với nội dung hiển thị trên từng slide, chia khối ~25 từ, đã kiểm tra biên dịch sạch (không lỗi LaTeX).

Compile cả 3 file đều không còn lỗi (`pdflatex` sạch, chỉ còn warning "Overfull \hbox" vô hại do font Courier trong `transcript.tex`).

## Biên soạn (compile)
```
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```
Chạy `pdflatex` 2 lần cuối để mục lục và trích dẫn (`\cite`) cập nhật đúng số trang/số thứ tự.

Với `slide.tex` và `transcript.tex`, chạy `pdflatex` tương tự (không cần `bibtex` với `transcript.tex` vì không có trích dẫn).
