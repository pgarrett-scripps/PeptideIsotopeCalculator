.PHONY: help venv install run-single run-multi run clean test lint format check

help:
	@echo "Available commands:"
	@echo "  make venv        - Create virtual environment with uv"
	@echo "  make install     - Install dependencies"
	@echo "  make run-single  - Run single isotope calculator app"
	@echo "  make run-multi   - Run multi isotope calculator app"
	@echo "  make run         - Run single isotope calculator app (default)"
	@echo "  make clean       - Remove cache and temporary files"
	@echo "  make lint        - Run code linting"
	@echo "  make format      - Format code with black"
	@echo "  make check       - Run linting and format check"

venv:
	uv venv

install:
	uv pip install -r requirements.txt

run-single:
	.venv/bin/streamlit run single_iso_app.py

run-multi:
	.venv/bin/streamlit run multi_iso_app.py

run: run-single

clean:
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".pytest_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".streamlit" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
	find . -type f -name "*~" -delete

lint:
	@command -v flake8 >/dev/null 2>&1 && .venv/bin/flake8 *.py || echo "flake8 not installed, skipping..."

format:
	@command -v black >/dev/null 2>&1 && .venv/bin/black *.py || echo "black not installed, skipping..."

check: lint
	@command -v black >/dev/null 2>&1 && .venv/bin/black --check *.py || echo "black not installed, skipping..."
