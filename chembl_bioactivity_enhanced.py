#!/usr/bin/env python3
# ===========================
# Neat ChEMBL Bioactivity Report (Integrated, Binder-hardened)
# ===========================
# Features
# - Pharmacodynamics from ChEMBL (human targets)
# - Search/sort/filter table (itables/DataTables)
# - Download CSV/Excel (separator selectable)
# - 2D structure (RDKit if available; else PubChem PNG)
# - 3D interactive structure (py3Dmol; requires PubChem 3D)
# - PubChem properties (robust walker for Experimental/Computed)
# - Binder/Voila-hardened networking (UA + retries + TXT fallback)
#
# Notes
# * RDKit via conda-forge is recommended (see environment.yml).
# * 3D model depends on PubChem having a 3D conformer.
# * Experimental properties are heterogeneous; availability varies.
# * This cell only builds the UI; no network calls occur until you press "Search".

import io, base64, json, re, requests, contextlib, warnings, time, unicodedata
from urllib.parse import quote
import pandas as pd
from IPython.display import display, clear_output, Markdown, HTML, Image
import ipywidgets as widgets

# Silence known deprecation warning from chembl_webresource_client
warnings.filterwarnings(
    "ignore",
    message="pkg_resources is deprecated as an API",
    category=UserWarning,
    module=r"chembl_webresource_client.*"
)

from chembl_webresource_client.new_client import new_client
from itables import init_notebook_mode, show, options as itbl_options
from chembl_bioactivity_enhanced import clean_activity_df, filter_by_threshold, pk_profile_for_compound

itbl_options.warn_on_undocumented_option = False
itbl_options.pageLength = 50
itbl_options.lengthMenu = [[10, 25, 50, 100, -1], [10, 25, 50, 100, "All"]]
init_notebook_mode(all_interactive=True)

# ----------------------------
# Interactive UI
# ----------------------------
def interactive_mode():
    text = widgets.Text(
        value='scopolamine',
        description='Compound:',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='420px')
    )
    act_filter = widgets.SelectMultiple(
        options=['IC50', 'Ki', 'KA', 'Kd'],
        value=['IC50', 'Ki'],
        description='Activity filter',
        layout=widgets.Layout(width='220px', height='110px')
    )
    sort_col = widgets.Dropdown(
        options=['Target', 'Activity', 'Value', 'Units', 'Kd (nM) (from KA)'],
        value='Target',
        description='Sort by:',
        layout=widgets.Layout(width='220px')
    )
    sort_asc = widgets.ToggleButtons(
        options=[('Asc', True), ('Desc', False)],
        value=True,
        description='Order:',
        layout=widgets.Layout(width='200px')
    )
    _sep_options = [('Comma ,', ','), ('Semicolon ;', ';'), ('Tab \\t', '\\t')]
    sep_choice = widgets.Dropdown(
        options=_sep_options,
        value=_sep_options[0][1],
        description='CSV separator:',
        layout=widgets.Layout(width='240px')
    )
    run_btn = widgets.Button(description="Search", button_style='primary')
    out = widgets.Output()

    controls = widgets.HBox([text, run_btn])
    filters  = widgets.HBox([act_filter, sort_col, sort_asc, sep_choice])

    def on_click(_):
        with out:
            clear_output()

            compound = text.value
            compound = unicodedata.normalize("NFKC", compound).strip()
            if not compound:
                print("🛑 Please enter a compound name.")
                return

            display(Markdown(f"## Results for **{compound}**"))
            display(Markdown("Data sources: **ChEMBL** (bioactivity), **PubChem** (structure & properties)."))

            # --- ChEMBL PD table
            df_pd = pd.DataFrame()
            try:
                chembl_id = new_client.molecule.filter(pref_name__iexact=compound)[0]['molecule_chembl_id']
                acts = new_client.activity.filter(
                    molecule_chembl_id=chembl_id,
                    target_organism__iexact='Homo sapiens'
                ).only([
                    'target_chembl_id',
                    'standard_type',
                    'standard_value',
                    'standard_units'
                ])
                df_pd = clean_activity_df(acts)

                if not df_pd.empty:
                    selected = list(act_filter.value)
                    if selected:
                        df_pd = df_pd[df_pd['Activity'].isin(selected)]
                    df_pd = filter_by_threshold(df_pd, 'Value_nM', max_nM=10000)
                    display(Markdown("### Pharmacodynamic Bioactivities (Homo sapiens)"))
                    show(df_pd, classes="display compact cell-border", maxBytes=0, pageLength=50)
                    display(HTML(make_download_link(df_pd, "bioactivity.csv", "csv", sep=sep_choice.value)))
                    display(HTML(make_download_link(df_pd, "bioactivity.xlsx", "xlsx")))
                else:
                    display(Markdown("> No human bioactivity rows returned by ChEMBL."))
            except Exception as e:
                display(Markdown(f"**ChEMBL error:** {e}"))

            # --- PubChem + Pharmacokinetics
            cid = None
            try:
                cid = new_client.molecule.filter(pref_name__iexact=compound)[0]['molecule_chembl_id']
            except:
                pass
            if cid:
                display(Markdown(f"### PubChem\nCID: **{cid}**"))
                df_pk = pk_profile_for_compound(cid)
                show(df_pk, classes="display compact cell-border", maxBytes=0, pageLength=50)
                display(HTML(make_download_link(df_pk, "pharmacokinetic_data.csv", "csv", sep=sep_choice.value)))
            else:
                # Pharmacokinetics SKIPPED
                pass

    run_btn.on_click(on_click)
    display(widgets.VBox([controls, filters, out]))

# Launch the UI
interactive_mode()
