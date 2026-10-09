# Worked seismology notebooks

Read the [rendered tutorials](https://denolle-lab.github.io/gesc-python/), or run the notebooks in this directory. They use synthetic observations at field scales; the Cook Inlet tutorial separately cites published KKFL-S acquisition facts. No continuous DAS data or field velocity model is bundled.

1. `01_critical_zone.ipynb`: regional depth sensitivity, a shallow softening event, correlated frequency errors, band-limited resolution, and synthetic time histories.
2. `02_volcano.ipynb`: shallow/deeper elastic-change hypotheses, posterior correlations, averaging matrices, and frequency coverage.
3. `03_scholte_water.ipynb`: analytical finite-water verification and water–solid nuisance tradeoffs.
4. `04_time_lapse_tomography.ipynb`: four-cell straight-ray acquisition and joint/two-stage inference with full covariance.
5. `05_cook_inlet_das.ipynb`: KKFL-S acquisition context, ideal DAS direction/gauge response, a constrained absolute-dispersion Vs fit, and seabed time-lapse inference with water uncertainty.

## Run locally

Use Python 3.12 for the pinned tutorial environment. From the checkout:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements/tutorials.txt
python -m pip install -e '.[dev,tutorials]'
python scripts/build_tutorials.py --output build/site
python -m http.server 8000 --bind 127.0.0.1 --directory build/site
```

Open `http://localhost:8000`. To work interactively, open a notebook in Jupyter or a notebook editor using that environment; start in the repository root or `tutorials/`. The downloadable executed notebooks also need a checkout (including the small `tutorials/` helpers), rather than only a standalone ipynb.

Allow a few minutes per notebook on a laptop: depth linearization re-solves several perturbed models with the dense reference backend.

The build runs each notebook from a fresh kernel, clears stored outputs, stops on any exception/assertion failure, checks numerical metrics and embedded figures, and checks local links. The optional `--write-executed` flag refreshes the committed notebook outputs. Pages publishes the freshly executed HTML, notebooks, figures, and `validation.json`. Dependency installation requires network access; the actual scientific calculations and tests are offline.

## Scope and interpretation

Depth coefficients perturb log Vs in each physical region, with Vp/density fixed. The halfspace is one coefficient over all its sensitivity below the last interface; it is not a resolved deep layer. Gaussian monitoring intervals condition on background structure, priors, covariances, and the direct fundamental phase observable. They do not convert coda MWCS/stretching automatically, marginalize background geology, or identify a unique physical process.

Tomography and DAS response helpers are teaching examples separate from the GESC core. The absolute Cook Inlet fit is a small MAP experiment with a local covariance approximation and fixed thickness/Vp/density. It does not implement a general structural inversion or reproduce a field Scholte profile. Explicit unknowns and source distinctions appear in each notebook.
