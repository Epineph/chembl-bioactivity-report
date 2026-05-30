# Chembl Bioactivity Report

## Run it in your browser

Click the **Launch Binder** badge below to run the notebook online in your
browser and investigate the pharmacodynamic activity of a molecule of interest.

The data sources are official and reliable. They are listed in both the rendered
notebook and the source code available in this repository. You do **not** need
to download or install anything to use the notebook.

Results are displayed in an interactive table and can also be downloaded as
Excel files (`.xlsx`) or delimited text files (`.csv`, semicolon-separated, or
tab-separated).

<p align="center">
  <strong>Launch Binder</strong><br>
  ⬇️<br>

</p>

<p align="center">
  <a href="https://mybinder.org/v2/gh/Epineph/chembl-bioactivity-report/main?urlpath=voila/render/chembl.ipynb">
    <img src="https://mybinder.org/badge_logo.svg" alt="Launch Binder">
  </a>
</p>
<!-- [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Epineph/chembl-bioactivity-report/main?urlpath=voila/render/chembl.ipynb)-->


**Neatly fetch human bioactivity (IC₅₀, Kᵢ, Kₐ, …) from ChEMBL and render as a GitHub‐style table.**\
[![Python ≥3.9](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)\
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

---

## 🚀 Features

- **Lookup** by compound name (case-insensitive)
- **Fetch** all _Homo sapiens_ bioactivities
- **Compute** Kd (nM) for association constants
- **Render** as Markdown table (CLI) or HTML (Jupyter)

## ⚙️ Installation

### Create a virtual environment

You can use any name for your virtual environment, below chembl-bioactivity is used as an example.

```bash
# Create a fresh environment dedicated to this workflow
micromamba create -n chembl-bioactivity -c conda-forge \
  python=3.11 \
  "numpy<2" \
  rdkit \
  chembl-webresource-client pandas tabulate ipywidgets itables openpyxl \
  requests pubchempy py3Dmol pillow -y

# Register a Jupyter kernel for the new env
micromamba activate chembl-bioactivity
python -m ipykernel install --user --name chembl-bioactivity --display-name "Python (chembl-bioactivity)"

```

### Install some dependencies

```bash
python -m pip install rdkit numpy pandas jupyterlab ipykernel ipywidgets itables pillow openpyxl voila
# If using mamba or conda rather than micromamba, use that instead
# Alternatively, use pip to install all the packages
micromamba install -c conda-forge chembl-webresource-client pubchempy py3Dmol tabulate

```

Some users may have to include _python3_ instead of _python_ in the command, so:

```bash
python3 -m pip install rdkit numpy pandas jupyterlab ipykernel ipywidgets itables pillow openpyxl voila chembl-webresource-client pubchempy py3Dmol tabulate
```

Then just launch jupyter lab (or use the binder badge above) locally:

```bash
jupyter lab
```

And open the _chembl_bioactivity.ipynb_
