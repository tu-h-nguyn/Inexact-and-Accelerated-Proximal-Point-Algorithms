# Inexact and Accelerated Proximal Point Algorithms

[![CI](https://github.com/tu-h-nguyn/inexact-and-accelerated-proximal-point-algorithms/actions/workflows/ci.yml/badge.svg)](https://github.com/tu-h-nguyn/inexact-and-accelerated-proximal-point-algorithms/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Tests](https://img.shields.io/badge/tests-10%20passing-brightgreen.svg)](tests/)
[![Bounds](https://img.shields.io/badge/bounds%20checked-38%20items-brightgreen.svg)](matlab/verify_bounds.m)
[![Octave](https://img.shields.io/badge/MATLAB%20%2F%20Octave-portable-orange.svg)](matlab/)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)

The proximal point method assumes you can solve

$$
\mathrm{prox}_{\lambda F}(y) = \arg\min_{z}\ F(z) + \tfrac{1}{2\lambda}\|z-y\|^2
$$

exactly. In practice you never can — the subproblem is itself solved by an
iterative method, stopped after finitely many steps. This report asks the
question that follows: **how inexactly can you solve it and still keep
Nesterov acceleration?**

The answer turns on *which* notion of inexactness you use. Two are analysed, and
they behave very differently: with the same error schedule, the accumulated error
term decays like $k^{-1/2}$ under the first and $k^{-2}$ under the second. Only
the second is fast enough to preserve the $O(1/k^2)$ rate.

**The bounds here are checked, not illustrated.** The central convergence bound is
tested at every outer iteration, for every inner-iteration schedule, and CI fails
the build if it is ever violated.

---

## What is in this repository

| | |
|---|---|
| **`matlab/`** | The LASSO experiment, portable between MATLAB and Octave. `verify_bounds.m` is the headless gate: 38 checks of the bounds, the theorem's own hypotheses, and the closed-form identities. |
| **`test_01/`** | A second, sharper experiment in Python on a 1-D quartic, where the predicted decay orders can be measured directly. Reproduces byte for byte. |
| **`tests/`** | 10 tests asserting the report's claims against the committed results. |
| **`content/`, `main.tex`** | The report (Vietnamese). `slide.tex` and `transcript.tex` reuse the same sources — see *One source, three outputs* below. |

## The two error criteria

Write $\hat x \approx \mathrm{prox}_{\lambda F}(y)$ with accuracy $\varepsilon$. The
two criteria differ in how the per-step error enters the accumulated term
$\delta_k$:

| | accumulation rule | measured decay |
|---|---|---|
| **Type 1** | $\eta_k = \sum_i \varepsilon_i/\alpha_i$ accumulates *before* being squared | $k^{-0{,}52}$ |
| **Type 2** | $\delta_{k+1} = (1-\alpha_k)\delta_k + \varepsilon_k^2/(2\lambda_k)$ — squared immediately | $k^{-1{,}99}$ |

Both measured at $q = 7/4$ over iterations 500–2000, against the predicted
$k^{-1/2}$ and $k^{-2}$. The type-2 slope lands within 0,6% of $-2$.

It would be tidy to say type 2 is simply the better bound. It is not, and the
code checks the real statement: at the first outer step the ratio
$\delta_2/\delta_1$ is **exactly** $(\lambda\rho+1)/(\lambda\rho) = 26$ — type 2 is
*26 times worse*. It crosses below type 1 only at step 9, and ends about 100×
tighter. "Recovering the rate" is an asymptotic statement, and the verification
script asserts the crossover rather than papering over it.

## What the verification actually checks

`make verify` runs the algorithm headless and checks 38 items:

- **The convergence bound of Theorem 3.2**, for all five inner-iteration
  schedules: **0 violations out of 121 outer steps**, each.
- **The theorem's own hypothesis.** It requires
  $a \le \alpha_k^2/\bigl((1-\alpha_k)A_k\lambda_k\bigr) \le 2$. The closed form used
  in `run_outer.m` makes that ratio **exactly 1** at every step — verified to
  $7\times10^{-16}$. A bound is only as good as the hypothesis it assumes, so the
  hypothesis is checked too.
- **Closed-form identities**, to machine precision:
  $\varepsilon_2/\varepsilon_1 = \sqrt{(\lambda\rho+1)/(\lambda\rho)}$, the recurrence
  $A_{k+1} = (1-\alpha_k)A_k$, the step-1 ratio above.
- **That inexactness matters.** Solving the subproblem with a single inner
  iteration leaves the final gap at $1{,}0\times10^{-1}$, almost seven orders of
  magnitude worse than $1{,}5\times10^{-8}$. A *growing* schedule ($T_k \sim \log k$ or
  $\sqrt{k}$) recovers essentially all of it; a coarse constant schedule does not.

The script is seed-independent — it passes on seeds 1, 7 and 42 — and it was
mutation-tested: perturbing the coefficient in the $\alpha_k$ formula from 4 to
3.5 is caught, and the build fails.

### Two experiments, two different standards of proof

The repository deliberately gates them differently, because they are not equally
reproducible:

- **`test_01/` (Python)** is a deterministic 1-D recursion — no RNG, no parallel
  reduction. It reproduces **byte for byte**, so CI does gate it with
  `git diff --exit-code`.
- **`matlab/`** cannot be gated that way. MATLAB and Octave ship different random
  number generators, so the same seed builds a *different problem instance*.
  Comparing digits would fail on the first environment change while catching no
  real error. `verify_bounds.m` instead checks the claims that hold for **every**
  instance, since those are what the theorems assert.

## Running it

```bash
make install     # Python dependencies
make test        # 10 tests
make verify      # 38 bound checks, headless Octave (no MATLAB licence needed)
make experiment  # re-run the Python experiment
make figures     # re-run the MATLAB/Octave experiment and redraw matlab/*.pdf
make report      # build main.pdf (needs a TeX distribution)
```

The MATLAB code runs unmodified under **GNU Octave**, so the numerical results
can be reproduced without a MATLAB licence. `export_fig_pdf.m` picks
`exportgraphics` or `print` depending on what the host provides.

## One source, three outputs

`main.tex` (report), `slide.tex` (beamer slides) and `transcript.tex` (speaking
script) all `\input` the same files from `content/`, and each shows a different
slice via the `comment` package:

| environment | report | slides | script |
|---|---|---|---|
| `essayonly` | shown | hidden | hidden |
| `slidesonly` | hidden | shown | hidden |
| `scriptonly` | hidden | hidden | shown |

Full proofs live in `essayonly`, so the slides keep only the statements while the
report keeps everything — and the two cannot drift apart, because there is one copy.

## Authors

Group project for HCMUS, supervised by TS. Nguyễn Đăng Khoa — Nguyễn Hoàng Tú,
Bùi Công Hoàng Vũ, Huỳnh Trung Kiên, Nguyễn Ngọc Diễm Quỳnh. See
[`CITATION.cff`](CITATION.cff).

## AI transparency

Claude (Anthropic) was used as an assistant on the engineering around the report,
not on its mathematics:

- **The report and its proofs are the group's.**
- **AI-assisted work is the tooling.** Splitting the duplicated experiment script,
  making the MATLAB code run under Octave, `verify_bounds.m`, the test suite,
  the packaging and CI, and this page.
- **Two real bugs were found and fixed in the process**, both noted in the commit
  history: `main_experiment.m` contained the entire script twice (the first copy
  was dead code, wiped by the second copy's `clear`), and `run_experiment.py`
  wrote its figures to a directory the LaTeX sources do not look in, so the paper
  could not find them.

## Licence

MIT for the code ([`LICENSE`](LICENSE)); the report text and figures are the
authors'.
