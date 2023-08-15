tidy:
	isort lig3dlens tests --profile black
	black lig3dlens tests

lint:
	ruff check lig3dlens tests
	isort lig3dlens tests --check-only --profile black
	black --diff lig3dlens tests

test:
	pytest -s --cov=lig3dlens --cov-report term-missing tests/
