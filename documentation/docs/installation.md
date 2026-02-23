# Installation

> Note: k4neo is still in development stage and can only be installed via a combination of conda, poetry and pip locally. We plan to provide k4neo on conda-forge and as container image

Clone the repository and change into the directory:


### Create a conda environment with all non-python dependencies.

```
conda env create -f k4neo.yaml -p k4neo_env
conda activate k4neo_env/
```

### Install k4neo package

```
poetry build -f wheel
pip install dist/k4neo-*-py3-none-any.whl
```

### 🧪 Run test suite

This will execute the comprehensive integration tests. You can use this to verify that your local installation works.

```
pytest --git-aware --symlink --stderr-bytes 100000 tests/
```