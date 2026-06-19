import io
import unittest

import pandas as pd

from chembl_bioactivity_enhanced import (
    active_metabolite_notes_for_compound,
    aggregate_activity_replicates,
    build_interactive_html_report,
    build_pdf_report,
    clean_activity_df,
    estimate_exposure_thresholds,
    filter_by_threshold,
    has_complete_pk_descriptors,
    limited_report_dataframe,
    merge_pk_descriptor_fallback,
    missing_pk_descriptor_names,
    pk_estimate_formula_items,
    pk_profile_from_descriptors,
    predict_pk_from_descriptors,
    report_row_limit,
    simulate_active_metabolite_curves,
    simulate_pk_curves,
)


class PharmacokineticPredictionTests(unittest.TestCase):
    def test_druglike_descriptors_produce_finite_pk_estimates(self):
        prediction = predict_pk_from_descriptors(
            {
                "MolecularWeight": 303.35,
                "XLogP": 1.2,
                "TPSA": 62.3,
                "HBondDonorCount": 1,
                "HBondAcceptorCount": 5,
                "RotatableBondCount": 4,
                "Charge": 0,
            }
        )

        self.assertGreater(prediction["fraction_absorbed_pct"], 50)
        self.assertLessEqual(prediction["oral_bioavailability_pct"], prediction["fraction_absorbed_pct"])
        self.assertGreater(prediction["vd_l_kg"], 0)
        self.assertGreater(prediction["clearance_ml_min_kg"], 0)
        self.assertGreater(prediction["half_life_h"], 0)
        self.assertIn("CNS penetration", prediction["bbb_class"])

    def test_polar_high_mw_compound_is_penalized(self):
        druglike = predict_pk_from_descriptors(
            {
                "MolecularWeight": 310,
                "XLogP": 2.2,
                "TPSA": 65,
                "HBondDonorCount": 1,
                "HBondAcceptorCount": 4,
                "RotatableBondCount": 4,
                "Charge": 0,
            }
        )
        polar = predict_pk_from_descriptors(
            {
                "MolecularWeight": 650,
                "XLogP": -1.0,
                "TPSA": 220,
                "HBondDonorCount": 8,
                "HBondAcceptorCount": 14,
                "RotatableBondCount": 16,
                "Charge": -1,
            }
        )

        self.assertGreater(druglike["fraction_absorbed_pct"], polar["fraction_absorbed_pct"])
        self.assertGreater(druglike["oral_bioavailability_pct"], polar["oral_bioavailability_pct"])
        self.assertGreater(druglike["logbb"], polar["logbb"])
        self.assertGreaterEqual(polar["lipinski_violations"], 3)

    def test_pk_profile_table_exposes_formula_context(self):
        df = pk_profile_from_descriptors(
            {
                "MolecularWeight": 310,
                "XLogP": 2.2,
                "TPSA": 65,
                "HBondDonorCount": 1,
                "HBondAcceptorCount": 4,
                "RotatableBondCount": 4,
                "Charge": 0,
            }
        )

        self.assertIn("Predicted terminal half-life", set(df["Parameter"]))
        self.assertIn("Evidence / formula", df.columns)
        self.assertIn("Confidence", df.columns)

    def test_profile_format_keeps_zero_and_trailing_zero_values(self):
        df = pk_profile_from_descriptors(
            {
                "MolecularWeight": 284.74,
                "XLogP": 3.0,
                "TPSA": 32.7,
                "HBondDonorCount": 0,
                "HBondAcceptorCount": 2,
                "RotatableBondCount": 1,
                "Charge": 0,
            }
        )

        hbd = df.loc[df["Parameter"] == "H-bond donors", "Estimate"].iloc[0]
        absorbed = df.loc[df["Parameter"] == "Predicted fraction absorbed orally", "Estimate"].iloc[0]

        self.assertEqual(hbd, "0")
        self.assertGreaterEqual(float(absorbed), 80)

    def test_descriptor_completeness_treats_zero_as_present(self):
        descriptors = {
            "MolecularWeight": 300,
            "XLogP": 0,
            "TPSA": 0,
            "HBondDonorCount": 0,
            "HBondAcceptorCount": 0,
            "RotatableBondCount": 0,
            "FormalCharge": 0,
        }

        self.assertTrue(has_complete_pk_descriptors(descriptors))
        self.assertFalse(has_complete_pk_descriptors({"MolecularWeight": 300, "XLogP": 0}))

    def test_pk_profile_does_not_predict_from_incomplete_descriptors(self):
        df = pk_profile_from_descriptors({"MolecularWeight": 300}, label="Partial")

        self.assertEqual(df.loc[0, "Parameter"], "Predicted pharmacokinetics unavailable")
        self.assertIn("Missing required descriptor", df.loc[0, "Evidence / formula"])
        self.assertIn("TPSA", missing_pk_descriptor_names({"MolecularWeight": 300}))

    def test_pk_prediction_rejects_incomplete_descriptors(self):
        with self.assertRaisesRegex(ValueError, "incomplete descriptors"):
            predict_pk_from_descriptors({"MolecularWeight": 300})

    def test_descriptor_fallback_fills_missing_values_without_overwriting(self):
        merged = merge_pk_descriptor_fallback(
            {
                "MolecularWeight": 300,
                "XLogP": 1.5,
                "DescriptorSource": "PubChem basic properties",
            },
            {
                "MolecularWeight": 310,
                "TPSA": 65,
                "HBondDonorCount": 1,
                "HBondAcceptorCount": 4,
                "RotatableBondCount": 3,
                "Charge": 0,
                "DescriptorSource": "RDKit from SMILES",
            },
        )

        self.assertEqual(merged["MolecularWeight"], 300)
        self.assertEqual(merged["XLogP"], 1.5)
        self.assertEqual(merged["TPSA"], 65)
        self.assertIn("RDKit from SMILES fallback", merged["DescriptorSource"])

    def test_simulation_returns_route_summary_and_curve(self):
        prediction = predict_pk_from_descriptors(
            {
                "MolecularWeight": 310,
                "XLogP": 2.2,
                "TPSA": 65,
                "HBondDonorCount": 1,
                "HBondAcceptorCount": 4,
                "RotatableBondCount": 4,
                "Charge": 0,
            }
        )
        curve, summary = simulate_pk_curves(
            prediction,
            dose_amount=10,
            dose_unit="milligram",
            routes=["Oral", "Intravenous"],
            mec_mg_l=0.01,
            mtc_mg_l=1.0,
            injection_concentration_mg_ml=5,
        )

        self.assertEqual(set(summary["Route"]), {"Oral", "Intravenous"})
        self.assertGreater(summary.loc[summary["Route"] == "Intravenous", "Cmax (mg/L)"].iloc[0], 0)
        self.assertEqual(summary.loc[summary["Route"] == "Intravenous", "Volume to administer (mL)"].iloc[0], 2)
        self.assertIn("Elimination t1/2 (h)", summary.columns)
        self.assertIn("Absorption t1/2 (h)", summary.columns)
        self.assertIn("Apparent terminal t1/2 (h)", summary.columns)
        self.assertIn("Terminal phase driver", summary.columns)
        self.assertIn("Post-peak 50% decline (h)", summary.columns)
        self.assertIn("Concentration (mg/L)", curve.columns)

    def test_route_half_times_separate_elimination_and_absorption(self):
        prediction = predict_pk_from_descriptors(
            {
                "MolecularWeight": 310,
                "XLogP": 2.2,
                "TPSA": 65,
                "HBondDonorCount": 1,
                "HBondAcceptorCount": 4,
                "RotatableBondCount": 4,
                "Charge": 0,
            }
        )
        _, summary = simulate_pk_curves(
            prediction,
            dose_amount=10,
            dose_unit="milligram",
            routes=["Oral", "Subcutaneous", "Intravenous"],
        )

        self.assertEqual(summary["Elimination t1/2 (h)"].nunique(), 1)
        self.assertNotEqual(
            summary.loc[summary["Route"] == "Oral", "Post-peak 50% decline (h)"].iloc[0],
            summary.loc[summary["Route"] == "Subcutaneous", "Post-peak 50% decline (h)"].iloc[0],
        )

    def test_apparent_terminal_half_life_can_be_absorption_limited_by_route(self):
        prediction = {
            "fraction_absorbed_pct": 90,
            "oral_bioavailability_pct": 75,
            "tpsa": 65,
            "rotb": 4,
            "mw": 310,
            "charge": 0,
            "vd_l_kg": 0.1,
            "clearance_ml_min_kg": 20,
        }

        _, summary = simulate_pk_curves(
            prediction,
            dose_amount=10,
            dose_unit="milligram",
            routes=["Oral", "Subcutaneous", "Intravenous"],
        )

        elimination_half_life = summary.loc[summary["Route"] == "Intravenous", "Elimination t1/2 (h)"].iloc[0]
        oral_terminal = summary.loc[summary["Route"] == "Oral", "Apparent terminal t1/2 (h)"].iloc[0]
        subcutaneous_terminal = summary.loc[
            summary["Route"] == "Subcutaneous", "Apparent terminal t1/2 (h)"
        ].iloc[0]

        self.assertGreater(oral_terminal, elimination_half_life)
        self.assertGreater(subcutaneous_terminal, oral_terminal)
        self.assertEqual(
            summary.loc[summary["Route"] == "Subcutaneous", "Terminal phase driver"].iloc[0],
            "Absorption-limited (flip-flop)",
        )

    def test_descriptor_formula_items_include_estimate_context(self):
        labels = {label for label, _ in pk_estimate_formula_items()}

        self.assertIn("Predicted fraction absorbed", labels)
        self.assertIn("Clark-style brain:blood estimate", labels)
        self.assertIn("Apparent terminal half-life", labels)
        self.assertIn("Post-peak half-time", labels)

    def test_reference_dose_can_estimate_mec_and_mtc(self):
        prediction = predict_pk_from_descriptors(
            {
                "MolecularWeight": 310,
                "XLogP": 2.2,
                "TPSA": 65,
                "HBondDonorCount": 1,
                "HBondAcceptorCount": 4,
                "RotatableBondCount": 4,
                "Charge": 0,
            }
        )

        mec, mtc, df = estimate_exposure_thresholds(
            prediction,
            active_dose_amount=5,
            toxic_dose_amount=50,
            reference_dose_unit="milligram",
            reference_route="Oral",
            cmax_fraction=0.5,
        )

        self.assertGreater(mec, 0)
        self.assertGreater(mtc, mec)
        self.assertEqual(set(df["Threshold"]), {"Estimated MEC", "Estimated MTC"})

    def test_tramadol_metabolite_curve_is_generated(self):
        prediction = predict_pk_from_descriptors(
            {
                "MolecularWeight": 263.38,
                "XLogP": 2.5,
                "TPSA": 32.7,
                "HBondDonorCount": 1,
                "HBondAcceptorCount": 2,
                "RotatableBondCount": 4,
                "Charge": 0,
            }
        )
        curve, summary = simulate_pk_curves(
            prediction,
            dose_amount=50,
            dose_unit="milligram",
            routes=["Oral"],
            body_weight_kg=70,
        )
        metabolite_curve, metabolite_summary = simulate_active_metabolite_curves(
            "tramadol", curve, summary, prediction, body_weight_kg=70
        )

        self.assertFalse(metabolite_curve.empty)
        self.assertIn("O-desmethyltramadol", metabolite_summary.loc[0, "Metabolite"])
        self.assertGreater(metabolite_summary.loc[0, "Relative potency vs parent"], 100)

    def test_active_metabolite_notes_include_codeine_caveat(self):
        df = active_metabolite_notes_for_compound("codeine")

        self.assertIn("Morphine", df.loc[0, "Active metabolite(s)"])
        self.assertIn("CYP2D6", df.loc[0, "Main pathway"])


class ActivityCleanupTests(unittest.TestCase):
    def test_activity_values_are_normalized_to_nm(self):
        df = clean_activity_df(
            [
                {
                    "target_chembl_id": "CHEMBL1",
                    "standard_type": "Ki",
                    "standard_value": "5",
                    "standard_units": "uM",
                },
                {
                    "target_chembl_id": "CHEMBL2",
                    "standard_type": "KA",
                    "standard_value": "1000000000",
                    "standard_units": "M^-1",
                },
            ]
        )

        self.assertEqual(df.loc[0, "Value_nM"], 5000)
        self.assertEqual(df.loc[1, "Value_nM"], 1)
        self.assertEqual(df.loc[1, "Kd (nM) (from KA)"], 1)

    def test_activity_cleanup_preserves_assay_context(self):
        df = clean_activity_df(
            [
                {
                    "target_pref_name": "Human receptor",
                    "target_chembl_id": "CHEMBL1",
                    "standard_type": "Ki",
                    "standard_relation": "<",
                    "standard_value": "5",
                    "standard_units": "nM",
                    "pchembl_value": "8.3",
                    "assay_chembl_id": "CHEMBLASSAY1",
                    "assay_type": "B",
                    "assay_description": "Binding assay",
                    "data_validity_comment": "Outside typical range",
                    "potential_duplicate": 1,
                }
            ]
        )

        self.assertEqual(df.loc[0, "Target"], "Human receptor")
        self.assertEqual(df.loc[0, "Relation"], "<")
        self.assertEqual(df.loc[0, "pChEMBL"], "8.3")
        self.assertEqual(df.loc[0, "Assay ChEMBL ID"], "CHEMBLASSAY1")
        self.assertEqual(df.loc[0, "Assay Type"], "B")
        self.assertEqual(df.loc[0, "Assay Description"], "Binding assay")
        self.assertEqual(df.loc[0, "Data Validity"], "Outside typical range")
        self.assertEqual(df.loc[0, "Potential Duplicate"], 1)

    def test_threshold_filter_preserves_unknown_units(self):
        df = clean_activity_df(
            [
                {
                    "target_chembl_id": "CHEMBL1",
                    "standard_type": "Ki",
                    "standard_value": "5",
                    "standard_units": "uM",
                },
                {
                    "target_chembl_id": "CHEMBL2",
                    "standard_type": "IC50",
                    "standard_value": "1",
                    "standard_units": "unknown",
                },
            ]
        )

        filtered = filter_by_threshold(df, max_nM=1000)

        self.assertEqual(list(filtered["Target"]), ["CHEMBL2"])

    def test_activity_replicates_are_averaged(self):
        df = clean_activity_df(
            [
                {
                    "target_chembl_id": "CHEMBL1",
                    "standard_type": "Ki",
                    "standard_value": "10",
                    "standard_units": "nM",
                },
                {
                    "target_chembl_id": "CHEMBL1",
                    "standard_type": "Ki",
                    "standard_value": "20",
                    "standard_units": "nM",
                },
            ]
        )

        aggregated = aggregate_activity_replicates(df)

        self.assertEqual(len(aggregated), 1)
        self.assertEqual(aggregated.loc[0, "Value_nM"], 15)
        self.assertEqual(aggregated.loc[0, "Units"], "nM")
        self.assertEqual(aggregated.loc[0, "Replicates"], 2)
        self.assertEqual(aggregated.loc[0, "Retained"], 2)
        self.assertEqual(aggregated.loc[0, "Excluded"], 0)
        self.assertAlmostEqual(aggregated.loc[0, "SD_nM"], 7.071, places=3)

    def test_activity_aggregation_keeps_assay_context_separate(self):
        df = clean_activity_df(
            [
                {
                    "target_chembl_id": "CHEMBL1",
                    "standard_type": "Ki",
                    "standard_relation": "=",
                    "standard_value": "10",
                    "standard_units": "nM",
                    "assay_chembl_id": "CHEMBLASSAY1",
                    "assay_type": "B",
                },
                {
                    "target_chembl_id": "CHEMBL1",
                    "standard_type": "Ki",
                    "standard_relation": "=",
                    "standard_value": "20",
                    "standard_units": "nM",
                    "assay_chembl_id": "CHEMBLASSAY1",
                    "assay_type": "B",
                },
                {
                    "target_chembl_id": "CHEMBL1",
                    "standard_type": "Ki",
                    "standard_relation": "=",
                    "standard_value": "100",
                    "standard_units": "nM",
                    "assay_chembl_id": "CHEMBLASSAY2",
                    "assay_type": "B",
                },
            ]
        )

        aggregated = aggregate_activity_replicates(df)
        assay_1 = aggregated[aggregated["Assay ChEMBL ID"] == "CHEMBLASSAY1"].iloc[0]
        assay_2 = aggregated[aggregated["Assay ChEMBL ID"] == "CHEMBLASSAY2"].iloc[0]

        self.assertEqual(len(aggregated), 2)
        self.assertEqual(assay_1["Value_nM"], 15)
        self.assertEqual(assay_1["Replicates"], 2)
        self.assertEqual(assay_2["Value_nM"], 100)
        self.assertEqual(assay_2["Replicates"], 1)

    def test_activity_replicate_mean_can_sigma_clip_outliers(self):
        df = clean_activity_df(
            [
                {
                    "target_chembl_id": "CHEMBL1",
                    "standard_type": "IC50",
                    "standard_value": value,
                    "standard_units": "nM",
                }
                for value in ["10", "11", "12", "100"]
            ]
        )

        aggregated = aggregate_activity_replicates(df, sigma_cutoff=1.0)

        self.assertEqual(aggregated.loc[0, "Value_nM"], 11)
        self.assertEqual(aggregated.loc[0, "Replicates"], 4)
        self.assertEqual(aggregated.loc[0, "Retained"], 3)
        self.assertEqual(aggregated.loc[0, "Excluded"], 1)

    def test_activity_aggregation_preserves_unknown_units(self):
        df = clean_activity_df(
            [
                {
                    "target_chembl_id": "CHEMBL1",
                    "standard_type": "Ki",
                    "standard_value": "5",
                    "standard_units": "unknown",
                },
                {
                    "target_chembl_id": "CHEMBL1",
                    "standard_type": "Ki",
                    "standard_value": "10",
                    "standard_units": "nM",
                },
            ]
        )

        aggregated = aggregate_activity_replicates(df)

        self.assertEqual(len(aggregated), 2)
        self.assertEqual(aggregated.loc[0, "Value_nM"], 10)
        self.assertTrue(pd.isna(aggregated.loc[1, "Value_nM"]))
        self.assertEqual(aggregated.loc[1, "Units"], "unknown")


class ReportExportTests(unittest.TestCase):
    def _tiny_png(self):
        from PIL import Image

        buf = io.BytesIO()
        Image.new("RGB", (4, 4), "white").save(buf, format="PNG")
        return buf.getvalue()

    def test_report_row_limit_normalizes_widget_values(self):
        self.assertEqual(report_row_limit(25), 25)
        self.assertEqual(report_row_limit("10"), 10)
        self.assertIsNone(report_row_limit("all"))
        self.assertIsNone(report_row_limit(0))
        self.assertEqual(report_row_limit(""), 50)
        self.assertEqual(report_row_limit("invalid"), 50)

    def test_limited_report_dataframe_respects_limit(self):
        df = pd.DataFrame({"Name": ["row-a", "row-b", "row-c"]})

        limited = limited_report_dataframe(df, row_limit=2)

        self.assertEqual(list(limited["Name"]), ["row-a", "row-b"])

    def test_html_report_limits_rows_and_escapes_content(self):
        df = pd.DataFrame({"Name": ["row-a", "row-b", "row-c"], "Value": ["<safe>", "B", "C"]})

        html_report = build_interactive_html_report(
            "Compound <Test>",
            [{"title": "Results <Table>", "df": df, "note": "Use <care>"}],
            row_limit=2,
        )

        self.assertIn("Compound &lt;Test&gt;", html_report)
        self.assertIn("Results &lt;Table&gt;", html_report)
        self.assertIn("Use &lt;care&gt;", html_report)
        self.assertIn("row-a", html_report)
        self.assertIn("row-b", html_report)
        self.assertNotIn("row-c", html_report)
        self.assertIn("(2 of 3 rows)", html_report)

    def test_html_report_embeds_image_sections(self):
        html_report = build_interactive_html_report(
            "Compound Test",
            [
                {
                    "title": "PK Curve",
                    "images": [{"title": "Concentration curve", "data": self._tiny_png(), "mime": "image/png"}],
                }
            ],
        )

        self.assertIn("PK Curve", html_report)
        self.assertIn("Concentration curve", html_report)
        self.assertIn("data:image/png;base64,", html_report)

    def test_pdf_report_builds_pdf_bytes_when_reportlab_is_available(self):
        try:
            import reportlab  # noqa: F401
        except Exception as exc:
            self.skipTest(f"reportlab unavailable: {exc}")

        df = pd.DataFrame({"Name": ["row-a", "row-b", "row-c"], "Value": [1, 2, 3]})

        pdf = build_pdf_report(
            "Compound Test",
            [
                {"title": "Results", "df": df},
                {"title": "PK Curve", "images": [{"title": "Concentration curve", "data": self._tiny_png()}]},
            ],
            row_limit=2,
        )

        self.assertTrue(pdf.startswith(b"%PDF-"))
        self.assertGreater(len(pdf), 1000)


if __name__ == "__main__":
    unittest.main()
