# Installation

## Local Development Checkout

Install the package in editable mode from the repository root:

```bash
python -m pip install -e .
```

The base install includes the database and the tools to query it, validate
records, and generate LaTeX tables, figures, and PowerPoint slides.

If you will edit code, stage records, or run tests and packaging checks,
install the contributor tools too:

```bash
python -m pip install -e ".[dev]"
```

For local documentation builds, install the documentation extra:

```bash
python -m pip install -e ".[docs]"
```

Both extras can be installed together:

```bash
python -m pip install -e ".[dev,docs]"
```

## Documentation Build

Build the documentation locally with:

```bash
python -m sphinx -b html docs docs/_build/html
```

Open `docs/_build/html/index.html` to read the built docs.

## Read the Docs

Read the Docs is configured through `.readthedocs.yaml`. The hosted build
installs the project with the `docs` optional dependency group, equivalent to:

```bash
python -m pip install ".[docs]"
```
