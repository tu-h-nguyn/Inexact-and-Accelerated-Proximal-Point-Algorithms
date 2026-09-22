PYTHON  ?= python3
OCTAVE  ?= octave --no-gui --quiet

.PHONY: help install test lint format verify experiment figures report slide transcript clean

help:
	@echo "install     install the Python dependencies for the numerical experiment"
	@echo "test        run the Python test suite"
	@echo "lint        run ruff over the Python sources"
	@echo "format      apply ruff's automatic fixes"
	@echo "verify      check the theoretical bounds numerically (Octave, headless)"
	@echo "experiment  re-run the Python experiment (rewrites test_01/*.csv and figures)"
	@echo "figures     re-run the MATLAB/Octave experiment and redraw matlab/*.pdf"
	@echo "report      build main.pdf with latexmk (needs a TeX distribution)"
	@echo "slide       build slide.pdf"
	@echo "transcript  build transcript.pdf"
	@echo "clean       remove build artifacts"

install:
	$(PYTHON) -m pip install -r test_01/requirements.txt pytest ruff

test:
	$(PYTHON) -m pytest

lint:
	$(PYTHON) -m ruff check .

format:
	$(PYTHON) -m ruff check --fix .

# The central claim of the report is a convergence bound. This runs the
# algorithm and checks that the bound is never violated, at every outer
# iteration, for every inner-iteration schedule.
verify:
	cd matlab && $(OCTAVE) --eval "exit(~verify_bounds())"

experiment:
	cd test_01 && $(PYTHON) run_experiment.py

figures:
	cd matlab && $(OCTAVE) --eval "main_experiment"

report:
	latexmk -pdf -interaction=nonstopmode main.tex

slide:
	latexmk -pdf -interaction=nonstopmode slide.tex

transcript:
	latexmk -pdf -interaction=nonstopmode transcript.tex

clean:
	rm -rf .pytest_cache .ruff_cache
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	latexmk -C || true
