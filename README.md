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

[![Python ≥3.9](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)\
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

---

## 🚀 Features

- **Lookup** by compound name (case-insensitive)
- **Fetch** all _Homo sapiens_ bioactivities
- **Compute** Kd (nM) for association constants
- **Estimate** structure-based pharmacokinetics for compounds without clinical PK data
- **Render** as Markdown table (CLI) or HTML (Jupyter)

## Predicted Pharmacokinetics

The app now adds a **Predicted Pharmacokinetics (structure-based)** table after
PubChem descriptor lookup. It estimates oral absorption potential, oral
bioavailability potential, BBB penetration tendency, volume of distribution,
clearance, and terminal half-life from molecular descriptors such as MW, logP,
TPSA, H-bond donors/acceptors, rotatable bonds, and formal charge.

The app also includes a dose/route simulator for exploratory concentration-time
curves. You can enter a dose in micrograms, milligrams, or grams; choose oral,
sublingual, intranasal, subcutaneous, and/or intravenous routes; enter optional
MEC/MTC thresholds in mg/L; and enter an injection concentration in mg/mL to
calculate `Volume to administer = Prescribed dose / Concentration` for IV or
subcutaneous routes. If MEC/MTC are unknown, you can provide a minimum active
dose and/or minimum toxic dose with a reference route; the app estimates the
threshold from a configurable fraction of the predicted reference-route Cmax.
The simulator reports Cmax, Tmax, AUC, half-life, onset to MEC, duration above
MEC, time above MTC, and route-specific assumptions.

The app adds a small curated active-metabolite/prodrug caveat table for common
examples such as codeine, diazepam, tramadol, clopidogrel, prednisone, and
enalapril. For selected cases where parent-only curves are especially misleading
(`tramadol`, `codeine`, and `diazepam`), it also plots a rough active-metabolite
or active-moiety curve. For tramadol, the M1 curve uses the parent/M1 MOR Ki
ratio as a potency-weighted approximation, so it should be interpreted as an
opioid-receptor active-moiety index rather than measured plasma concentration.

These are screening estimates for hypothesis generation, not clinical dosing or
safety guidance. The most evidence-aligned pieces are descriptor rules such as
Lipinski RO5, Veber oral exposure criteria, and the Clark-style logBB
relationship using logP and polar surface area. Clearance and half-life are
lower-confidence heuristics because accurate values usually require in vitro
metabolism/transport data or human/animal PK observations.

The concentration curves use a one-compartment model with first-order absorption
for non-IV routes and IV bolus kinetics for the intravenous route. They assume a
healthy young adult with normal renal and CYP activity and no major enzyme
inducers or inhibitors. These outputs are not dosing, prescribing, or safety
guidance. Active-metabolite models are curated heuristics only; the app does not
automatically infer enzyme pathways, stereochemistry, metabolite potency, or
genotype-specific metabolism for arbitrary molecules.

You can also call the estimator directly from Python:

```python
from chembl_bioactivity_enhanced import pk_profile_for_compound

# PubChem CID
df_pk = pk_profile_for_compound(3000322)

# Or a SMILES string, if RDKit is installed
df_pk = pk_profile_for_compound("CN1C(=O)N(C)c2ncn(C)c2C1=O")
```

## ⚙️ Installation

### Create a virtual environment

You can use any name for your virtual environment, below chembl-bioactivity is used as an example.

```bash
# Create a fresh environment dedicated to this workflow
micromamba create -n chembl-bioactivity -c conda-forge \
  python=3.11 \
  "numpy<2" \
  rdkit \
  chembl-webresource-client pandas tabulate ipywidgets itables openpyxl matplotlib \
  requests pubchempy py3Dmol pillow -y

# Register a Jupyter kernel for the new env
micromamba activate chembl-bioactivity
python -m ipykernel install --user --name chembl-bioactivity --display-name "Python (chembl-bioactivity)"

```

### Install some dependencies

```bash
python -m pip install rdkit numpy pandas matplotlib jupyterlab ipykernel ipywidgets itables pillow openpyxl voila
# If using mamba or conda rather than micromamba, use that instead
# Alternatively, use pip to install all the packages
micromamba install -c conda-forge chembl-webresource-client pubchempy py3Dmol tabulate

```

Some users may have to include _python3_ instead of _python_ in the command, so:

```bash
python3 -m pip install rdkit numpy pandas matplotlib jupyterlab ipykernel ipywidgets itables pillow openpyxl voila chembl-webresource-client pubchempy py3Dmol tabulate
```

Then just launch jupyter lab (or use the binder badge above) locally:

```bash
jupyter lab
```

And open the _chembl_bioactivity.ipynb_
