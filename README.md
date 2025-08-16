# Chembl Bioactivity Report

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Epineph/chembl-bioactivity-report/HEAD?urlpath=voila/render/chembl_bioactivity.ipynb)

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

```bash
pip install chembl-webresource-client pandas tabulate voila 
```
