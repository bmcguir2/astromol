# Example Notebooks And Google Colab

You can run these notebooks:

- locally, from a source checkout or installed package
- in Google Colab, opened directly from GitHub

The notebooks are maintained on GitHub. In Colab, the first setup cell installs
`astromol` from the `refactor` branch. The setup cells will need to switch to
PyPI when the package is released.

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

The links below open the `.ipynb` files in Colab or GitHub. Sphinx does not
render the notebooks, so building the docs does not require Pandoc.

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
