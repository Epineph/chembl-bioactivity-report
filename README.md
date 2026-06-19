# ChEMBL Bioactivity and Predicted Pharmacokinetics

Interactive Voila app for exploring molecule bioactivity, PubChem descriptors,
and transparent structure-based pharmacokinetic estimates.

The app is designed for hypothesis generation and teaching. It is **not** a
clinical dosing, prescribing, or safety tool.

## Launch Online

<p align="center">
  <a href="https://mybinder.org/v2/gh/Epineph/chembl-bioactivity-report/main?urlpath=voila/render/chembl.ipynb">
    <img src="https://mybinder.org/badge_logo.svg" alt="Launch Binder">
  </a>
</p>

[![Python >=3.11](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

## What It Does

- Looks up a compound by name.
- Fetches human bioactivity data from ChEMBL.
- Can average repeated target/activity bioactivity rows after nM normalization,
  with optional standard-deviation clipping.
- Resolves PubChem CID and computed descriptors.
- Estimates structure-based ADME/PK properties from descriptors.
- Simulates concentration-time curves by route and dose.
- Estimates MEC/MTC from optional reference active/toxic doses.
- Shows route-specific Cmax, Tmax, AUC, onset, duration, and post-peak decline.
- Provides curated active-metabolite caveats and limited active-moiety overlays.
- Exports tables and curves as CSV/XLSX, plus a PDF report and a collapsible standalone HTML report.

## Current Model Scope

The PK layer uses transparent heuristics rather than fitted clinical models.
Inputs include molecular weight, logP, TPSA, H-bond donors/acceptors, rotatable
bonds, and formal charge.

The strongest parts are descriptor rules and relationships such as Lipinski RO5,
Veber oral exposure criteria, and a Clark-style logBB estimate. Clearance,
volume of distribution, bioavailability, and half-life are lower-confidence
screening estimates because accurate values usually require in vitro metabolism,
transporter data, or observed animal/human PK.

The route simulator uses a one-compartment model with first-order absorption for
oral, sublingual, intranasal, and subcutaneous routes, and IV bolus kinetics for
intravenous route. It assumes a healthy 70 kg young adult unless changed in the
UI, normal renal/CYP activity, and no major enzyme inhibitors or inducers.

## MEC/MTC Calibration

If you know a rough minimum active dose or minimum toxic dose for a reference
route, the app can estimate concentration thresholds:

```text
estimated threshold = Ref Cmax fraction * predicted reference-dose Cmax
```

For example, `Ref Cmax fraction = 0.50` means MEC/MTC is estimated as half of
the predicted peak concentration for the supplied reference dose and route. This
is only a rough calibration tool; a real therapeutic window requires empirical
exposure-response data.

## Half-Life Interpretation

Elimination half-life is compound/system dependent and generally does not change
by route in the current one-compartment model. Route effects are shown through:

- absorption half-life
- apparent terminal half-life when absorption is slower than elimination
- Tmax
- onset and duration above MEC
- post-peak 50% decline time

Slow absorption can make the apparent post-peak decline longer than the intrinsic
elimination half-life. In flip-flop cases, the terminal slope is absorption
limited rather than elimination limited.

## Active Metabolites

The app includes curated caveats for selected examples such as codeine,
tramadol, diazepam, clopidogrel, prednisone, and enalapril. For a few compounds
where parent-only curves can be especially misleading, it also plots a rough
active-metabolite or active-moiety index.

These overlays are curated heuristics. The app does not automatically infer
enzyme pathways, stereochemistry, metabolite potency, or genotype-specific
metabolism for arbitrary molecules.

## Report Export

The app keeps individual CSV/XLSX downloads and also creates two summary reports:

- PDF report for printable table and static graph summaries.
- Standalone HTML report with clickable/collapsible sections and embedded static graphs.

Use the `Report rows/table` control to choose how many rows each displayed table
and report section should include. Dense concentration-time curves remain CSV
downloads rather than full report tables. The current report graphs are embedded
as PNG images, which work in both PDF and standalone HTML.

For pharmacodynamic bioactivity rows, `Mean repeats` averages normalized
`Value_nM` values within matching target/activity/assay/relation context.
This avoids averaging censored or assay-heterogeneous measurements together.
`Trim SD` removes values farther than the selected number of sample standard
deviations from that group mean before recomputing the displayed mean; `0`
disables clipping.

## Repository Layout

```text
chembl.ipynb                      Voila launcher notebook used by Binder
chembl_bioactivity_integrated.py  Interactive app and rendering logic
chembl_bioactivity_enhanced.py    Bioactivity cleanup, PK formulas, simulations
tests/                            Unit tests for PK helpers
binder/                           Binder environment and build hooks
requirements.txt                  pip dependencies for local installs
voila.json                        Voila configuration
```

## Local Installation

Conda/micromamba is recommended because RDKit is easiest to install from
conda-forge.

```bash
micromamba create -n chembl-pk -c conda-forge \
  python=3.11 rdkit numpy pandas matplotlib jupyterlab ipykernel ipywidgets \
  itables pillow openpyxl reportlab voila requests pubchempy py3Dmol tabulate -y

micromamba activate chembl-pk
python -m ipykernel install --user --name chembl-pk --display-name "Python (chembl-pk)"
```

Or install from `requirements.txt` in an existing environment:

```bash
python -m pip install -r requirements.txt
```

## Run Locally

Launch JupyterLab and open `chembl.ipynb`:

```bash
jupyter lab
```

Or render directly with Voila:

```bash
voila chembl.ipynb
```

## Run Tests

```bash
python -m unittest discover -s tests
```

## Data Sources

- ChEMBL: human target bioactivity data.
- PubChem: compound identifiers, structure, descriptors, images, and optional property scraping.
- RDKit: local descriptor calculation and 2D structure rendering when available.

## Important Disclaimer

All pharmacokinetic values are approximate screening estimates. They should not
be used to decide dose, route, medical treatment, toxicity risk, or personal drug
use. Clinical PK requires validated experimental data and professional review.
