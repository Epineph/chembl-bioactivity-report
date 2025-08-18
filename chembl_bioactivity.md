# %% [markdown]

# notebook to fetch pharmacodynamic bioactivity of compounds

Neat ChEMBL Bioactivity Report

- **Fetch pharmacodynamic bioactivity**
- Pulls all human bioactivities (*IC50*, *Ki*, *KA*, etc.) for a given compound.
- Looks up each target’s *preferred name*.
- **Builds a DataFrame and prints it as a pretty GitHub‑style Markdown table** when run
  in a terminal, or as an HTML table in Jupyter.
- When run from terminal, if using linux, make sure the script is executable by: sudo
  chmod +x chembl_bioactivity.py anr run from terminal:
- python chembl_bioactivity.py <compound>

______________________________________________________________________

# %% [markdown]

## 2) Install the required libraries

The notebook **requires the following libraries:**

- chembl-webresource-client
- pandas
- tabulate
- voila
- ipywidgets
- itables
- openpyxl

# Running script locally with jupyterlab notebook

Installing required packages for running it through *jupyter lab* locally using a localhost httml link, or host it on your own server with *jupyterhub*. 

If you run it locally on your pc, remember to get get the packages:
- *jupyterlab*
- *jupyter*
- *notebook* 
- *ipython*
- *ipykernel* 

This list above **is not exhaustive.**

It is recommended to install the packages within a virtual environment, using *venv* or a **cross-platform, language-agnostic binary package manager** like *conda*, *mamba* or *micromamba*:


# Installing packages
- create and activate the virtual environment you want to use
- make sure the python version is compatible

For example, you may create a virtual environment using *Python 3.11* with micromamba:

```bash
micromamba create -n testEnv -c conda-forge python=3.11
micromamba activate testEnv
```

# Confirm you are using your activated virtual environment

## Confirm version of python version in your virtual environment

After you have activated your environment, to confirm you are using it, you can use commands like:

```bash
python --version
```

The output from my test shows:

(testEnv)
~ [🐍 v3.11.13][🅒 testEnv]
❯ python --version
Python 3.11.13

So, after *confirming the python version matches the one used in your virtual environment*, **check the python location on disc that called from the terminal command when running the following command**:

```bash
which python
```

Which should show you are running python in a folder such as mamba or conda within your home folder, from my test, this is the output:


~ [🐍 v3.11.13][🅒 testEnv]
❯ which python
$HOME/.local/share/mamba/envs/testEnv/bin/python
(testEnv)

The path will probably be shown in the output as:
- */home/<your_user_name>/.local/share/mamba/envs/testEnv/bin/pyton*

This is equivalent to both:
- */home/$USER/.local/share/mamba/envs/testEnv/bin/pyton*
- *$HOME/.local/share/mamba/envs/testEnv/bin/pyton*

In fact, this is easy to confirm:

```bash
#/usr/bin/env bash
s1="/home/$USER/.local/share/mamba/envs/testEnv/bin/python"
s2="$HOME/.local/share/mamba/envs/testEnv/bin/python"

if [ "$(realpath -m "$s1")" = "$(realpath -m "$s2")" ]; then
    echo "Paths are identical: $s1 and $s2"
else
    echo "Paths differ: $s1 vs $s2"
fi
```

The output in this case would start with the home folder followed by your username, which I will just use **~/** below, but signifies the same:

**Paths are identical: ~/.local/share/mamba/envs/testEnv/bin/python and ~/.local/share/mamba/envs/testEnv/bin/python**

Of course, this is becaue */home/$USER"* and *$HOME* is the same thing. 

```bash
#/usr/bin/env bash
s3="/home/$USER"
s4="$HOME"

if [ "$(realpath -m "$s3")" = "$(realpath -m "$s4")" ]; then
    echo "Paths are identical: $s3 and $s4"
else
    echo "Paths differ: $s3 vs $s4"
fi

# This would echo "Paths differ: $3 vs $4"

# If s3="$HOME/$USER" and s4="$HOME"
# Then the output would echo that paths are unidentical.
# e.g.,:

s5="$HOME/$USER/.local/share/mamba/envs/testEnv/bin/python"
s6="$HOME/.local/share/mamba/envs/testEnv/bin/python"

if [ "$(realpath -m "$s5")" = "$(realpath -m "$s6")" ]; then
    echo "Paths are identical: $s5 and $s6"
else
    echo "Paths differ: $s5 vs $s6"
fi
```

The important thing, within the virtual environment, running **which python** (or with windows powershell: **where.exe python**) would point to a python in a subdirectory to mamba, conda or micromamba, which most often are subdirectories of your $HOME folder.

## Getting essential packages

If running the jupyter file locally from jupyterlab, make sure you have some essential packages like: *ipykernel*, *jupyterlab*, *jupyter*, *notebook* and *ipython*.

Install, making sure you are within the virtual environment, and using **python -m pip install <packages1> <package2> ... <package_n>, where *n = the number of the package appearing last* in the listed packages to be installed**.

## Packages crucial to run the script

Install with pip, preferably in your virtual environment, e.g.,:

```bash
# 
python -m pip install chembl-webresource-client pandas tabulate voila ipywidgets itables openpyxl
```

Or use your favorite packpage manager, e.g., **conda, mamba or micromamba**:

```bash
conda install -c conda-forge chembl-webresource-client pandas tabulate voila ipywidgets itables openpyxl
```

or clone the repository locally, cd to it, and run in your terminal (within your virtual env):

```bash
python -m pip install -r requirements.txt
```
# Script structure

## Imports & Setup

```python
import pandas as pd
from chembl_webresource_client.new_client import new_client
from IPython.display import display, clear_output, Markdown, HTML
import ipywidgets as widgets
from itables import init_notebook_mode, show
import io, base64

init_notebook_mode(all_interactive=True)
```

## Core helper functions

To find the **chembl_id** for the **compound name** string:

```python
# --- Core helper functions ---
def get_chembl_id(compound: str) -> str:
    mol_client = new_client.molecule
    res = mol_client.filter(pref_name__iexact=compound)
    if not res:
        raise ValueError(f"No ChEMBL entry for '{compound}'")
    return res[0]['molecule_chembl_id']
```

Fetch the activities we want to look at: 

```python
def fetch_activities(chembl_id: str) -> list[dict]:
    """Fetch all Homo sapiens bioactivities for the given ChEMBL ID."""
    act_client = new_client.activity
    acts = act_client.filter(
        molecule_chembl_id=chembl_id,
        target_organism__iexact='Homo sapiens'
    ).only([
        'target_chembl_id',
        'standard_type',
        'standard_value',
        'standard_units'
    ])
    return list(acts)
```

Mapping the **chembl_id** to the *preferred name* of the compound:

```python
def fetch_target_names(target_ids: set[str]) -> dict[str, str]:
    """Map each ChEMBL target ID to its preferred name."""
    tgt_client = new_client.target
    names = {}
    for tid in target_ids:
        rec = tgt_client.filter(target_chembl_id=tid).only(['pref_name'])
        names[tid] = rec[0]['pref_name'] if rec else tid
    return names
```

## Creating the dataframe and compute pharmacodynamic parameters

- **IQ50** – Concentration producing 50% of the maximum effect.  
- **Ki** – Inhibition constant; lower values indicate stronger binding to the target.  
- **KA** – Association (binding) constant; higher values indicate stronger affinity.  
- **Kd** – Dissociation constant; concentration at which half of the receptors are occupied.

```python
def build_activity_df(acts: list[dict]) -> pd.DataFrame:
    """Build a tidy DataFrame and compute Kd (nM) for KA entries."""
    rows = []
    for a in acts:
        tid   = a.get('target_chembl_id') or 'Unknown'
        typ   = a.get('standard_type')    or ''
        val   = a.get('standard_value')   or ''
        unit  = a.get('standard_units')   or ''
        kd_nm = ''
        # Compute Kd in nanomolar for association constants
        if typ.upper() == 'KA' and val and unit.strip() in ['M^-1','M-1','1/M']:
            try:
                kd_m  = 1.0 / float(val)
                kd_nm = round(kd_m * 1e9, 2)
            except Exception:
                kd_nm = ''
        rows.append({
            'Target (ChEMBL)': tid,
            'Activity':        typ,
            'Value':           val,
            'Units':           unit,
            'Kd (nM)':         kd_nm
        })
    df = pd.DataFrame(rows)
    # Replace ChEMBL IDs with human-readable names
    unique_tids = set(df['Target (ChEMBL)'])
    name_map = fetch_target_names(unique_tids)
    df['Target (ChEMBL)'] = df['Target (ChEMBL)'].map(name_map)
    # Keep only non-empty bioactivity entries
    df = df[['Target (ChEMBL)', 'Activity', 'Value', 'Units', 'Kd (nM)']]
    df.dropna(how='all', subset=['Value'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df
```

## Export table from compound search and download

Implementing a download link so that the table can be exported and downloaded as **.csv** or **.xlsx** (excel worksheet):

```python
# --- Utility: create browser download link ---
def make_download_link(df, filetype="csv", sep=","):
    if filetype == "csv":
        buf = io.StringIO()
        df.to_csv(buf, index=False, sep=sep)
        data = buf.getvalue()
        mime = "text/csv"
        filename = "results.csv"
    elif filetype == "xlsx":
        buf = io.BytesIO()
        df.to_excel(buf, index=False)
        data = buf.getvalue()
        mime = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        filename = "results.xlsx"
    else:
        raise ValueError("Unsupported type")

    b64 = base64.b64encode(
        data.encode() if isinstance(data, str) else data
    ).decode()
    return f'<a download="{filename}" href="data:{mime};base64,{b64}">⬇️ Download {filename}</a>'
```

## Widgets

Implementing interactive widgets, so that the user can interact with script in an intuitive way. That includes being able to write the name of the compound in the textbox, filter search and export and download it to your pc.

```python
# --- Interactive widget mode ---
def interactive_mode():
    text = widgets.Text(
        value='scopolamine',
        description='Compound:',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='400px')
    )
    button = widgets.Button(description="Search")
    filter_box = widgets.SelectMultiple(
        options=['IC50', 'Ki', 'KA', 'Kd'],
        value=['IC50', 'Ki'],
        description='Activity filter',
        layout=widgets.Layout(width='200px')
    )
    output = widgets.Output()

    def on_click(b):
        with output:
            clear_output()
            try:
                compound = text.value
                chembl_id = get_chembl_id(compound)
                acts = fetch_activities(chembl_id)
                df = build_activity_df(acts)
                selected = list(filter_box.value)
                if selected:
                    df = df[df['Activity'].isin(selected)]
                display(Markdown(f"""
                ### Results for **{compound}**  
                Data retrieved from [ChEMBL](https://www.ebi.ac.uk/chembl/).  
                Values aggregated per target/activity.
                """))
                show(df, classes="display compact cell-border", maxBytes=0)
                # Show safe browser download links
                display(HTML(make_download_link(df, "csv", sep=",")))
                display(HTML(make_download_link(df, "xlsx")))
            except Exception as e:
                print(f"❌ Error: {e}")

    button.on_click(on_click)

    display(widgets.VBox([
        text,
        button,
        filter_box,
        output
    ]))

interactive_mode()
```

# 🧪 Neat ChEMBL Bioactivity Report

- Pulls all human bioactivities (IC50, Ki, KA, etc.) for a given compound.
- Looks up each target’s preferred name.
- Builds a DataFrame and prints it as a pretty GitHub-style Markdown table\
  when run in a terminal, or as an HTML table in Jupyter.
