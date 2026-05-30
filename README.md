# Chembl Bioactivity Report

## Run it in your browser

Press on the icon below - **_Launch Binder_** to run the script online in your browser and investigate the pharmacodynamic action of the molecule of your interest. Sources, official and reliable, are listed both in the rendered script and visible in the raw code on github. You do not need to **download anything to use the script!**. 

Results are, in addition to being shown as you search for a given molecule in a table, available to download as *Excel* files, i.e., **.xlsx** format or *comma-seperated values*, i.e., **.csv** file (you can also choose to get the file as *semi-colon* (**;**) seperated values or *tab* (**/t**) seperated values if needed if you need it for a script that uses that).


[![Binder](https://mybinder.org/badge_logo.svg)](
  https://mybinder.org/v2/gh/Epineph/chembl-bioactivity-report/main?urlpath=voila/render/chembl.ipynb
)


<!-- H1 title; two spaces at end ⇒ hard line break for subtitle -->

**Neatly fetch human bioactivity (IC₅₀, Kᵢ, Kₐ, …) from ChEMBL and render as a
GitHub‐style table.**\
[![Python ≥3.9](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)\
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

______________________________________________________________________

## 🚀 Features

- **Lookup** by compound name (case-insensitive)
- **Fetch** all *Homo sapiens* bioactivities
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

Some users may have to include *python3* instead of *python* in the command, so:


```bash
python3 -m pip install rdkit numpy pandas jupyterlab ipykernel ipywidgets itables pillow openpyxl voila chembl-webresource-client pubchempy py3Dmol tabulate
```

Then just launch jupyter lab (or use the binder badge above) locally:

```bash
jupyter lab
```

And open the *chembl_bioactivity.ipynb* file from your browser.
