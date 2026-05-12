# Example Notebooks And Google Colab

The notebooks in this section are designed to work in two places:

- locally, from a source checkout or installed package
- in Google Colab, opened directly from GitHub

GitHub remains the source of truth. Colab is only the execution environment.
When a notebook is opened in Colab, its first setup cell installs `astromol`
from the current `refactor` branch. After the package is released on PyPI, those
setup cells should be changed to install `astromol` directly from PyPI.

## Opening A Notebook In Colab

Each notebook starts with an **Open in Colab** badge. The underlying URL follows
this pattern:

```text
https://colab.research.google.com/github/bmcguir2/astromol/blob/refactor/docs/notebooks/<notebook>.ipynb
```

In Colab, run the setup cell first, then run the remaining cells in order.
For the most common tasks, including the latest cumulative detections figure
and molecule slides, see [](faq.md).

## Notebook Index

The notebooks are stored as `.ipynb` files in GitHub and linked here rather
than rendered directly by Sphinx. This avoids requiring a Pandoc system
dependency in local and CI documentation builds.

| Notebook | Open in Colab | View on GitHub |
| --- | --- | --- |
| Quickstart | [Open in Colab](https://colab.research.google.com/github/bmcguir2/astromol/blob/refactor/docs/notebooks/01_quickstart.ipynb) | [01_quickstart.ipynb](https://github.com/bmcguir2/astromol/blob/refactor/docs/notebooks/01_quickstart.ipynb) |
| Reproduce figures | [Open in Colab](https://colab.research.google.com/github/bmcguir2/astromol/blob/refactor/docs/notebooks/02_reproduce_figures.ipynb) | [02_reproduce_figures.ipynb](https://github.com/bmcguir2/astromol/blob/refactor/docs/notebooks/02_reproduce_figures.ipynb) |
| Custom filtered views | [Open in Colab](https://colab.research.google.com/github/bmcguir2/astromol/blob/refactor/docs/notebooks/03_custom_views.ipynb) | [03_custom_views.ipynb](https://github.com/bmcguir2/astromol/blob/refactor/docs/notebooks/03_custom_views.ipynb) |
| Tables and slides | [Open in Colab](https://colab.research.google.com/github/bmcguir2/astromol/blob/refactor/docs/notebooks/04_tables_and_slides.ipynb) | [04_tables_and_slides.ipynb](https://github.com/bmcguir2/astromol/blob/refactor/docs/notebooks/04_tables_and_slides.ipynb) |

## Maintenance Notes

- Keep notebooks small and task-focused.
- Do not commit generated PDFs, PNGs, LaTeX fragments, or PowerPoint files.
- Prefer examples that call the same public functions documented in
  [](census-outputs.md).
- When the default branch changes from `refactor` to `main`, update the Colab
  badge URLs and setup cells.
- After PyPI release, update setup cells from GitHub install commands to
  `python -m pip install -q astromol`.
