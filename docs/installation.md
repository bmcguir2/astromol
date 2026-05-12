# Installation

## Local Development Checkout

Install the package in editable mode from the repository root:

```bash
python -m pip install -e .
```

The base package includes the database loader, census views, validation command,
curation staging workflow, LaTeX/table helpers, figure generation, and
PowerPoint slide generation.

For development and verification work, install the test runner too:

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

The generated HTML entry point is `docs/_build/html/index.html`.

## Read the Docs

Read the Docs is configured through `.readthedocs.yaml`. The hosted build
installs the project with the `docs` optional dependency group, equivalent to:

```bash
python -m pip install ".[docs]"
```
