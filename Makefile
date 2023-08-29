tidy:
	isort lig3dlens tests setup.py --profile black
	black lig3dlens tests setup.py

lint:
	ruff check lig3dlens tests setup.py
	isort lig3dlens tests setup.py --check-only --profile black
	black --diff lig3dlens tests setup.py

test:
	pytest -s --cov=lig3dlens --cov-report term-missing tests/
