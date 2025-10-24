.PHONY: help install install-dev test clean build upload upload-test docs

help: ## Show this help message
	@echo 'Usage: make [target]'
	@echo ''
	@echo 'Targets:'
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "  %-15s %s\n", $$1, $$2}' $(MAKEFILE_LIST)

install: ## Install the package
	pip install .

install-dev: ## Install the package in development mode
	pip install -e .[dev]

test: ## Run tests
	pytest

test-cov: ## Run tests with coverage
	pytest --cov=pymhm --cov-report=html --cov-report=term

clean: ## Clean build artifacts
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .pytest_cache/
	rm -rf .coverage
	rm -rf htmlcov/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

build: clean ## Build the package
	python setup.py sdist bdist_wheel

upload-test: build ## Upload to TestPyPI
	twine upload --repository testpypi dist/*

upload: build ## Upload to PyPI
	twine upload dist/*

docs: ## Build documentation
	cd docs && make html

lint: ## Run linting
	flake8 pymhm/
	black --check pymhm/
	mypy pymhm/

format: ## Format code
	black pymhm/

check: lint test ## Run all checks

release: clean build upload ## Build and upload to PyPI
