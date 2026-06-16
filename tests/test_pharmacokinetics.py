import unittest

from chembl_bioactivity_enhanced import (
    clean_activity_df,
    filter_by_threshold,
    pk_profile_from_descriptors,
    predict_pk_from_descriptors,
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
