create-env:
	conda env create -f environment.yml
	conda activate lig3dlens

tidy:
	isort . --profile black
	black .

test:
	pytest tests/