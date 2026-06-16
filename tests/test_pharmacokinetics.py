import unittest

from chembl_bioactivity_enhanced import (
    active_metabolite_notes_for_compound,
    clean_activity_df,
    estimate_exposure_thresholds,
    filter_by_threshold,
    pk_profile_from_descriptors,
    predict_pk_from_descriptors,
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
        self.assertIn("Concentration (mg/L)", curve.columns)

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


if __name__ == "__main__":
    unittest.main()
