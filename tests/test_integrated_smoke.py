import importlib
import sys
import types
import unittest
from unittest import mock

import pandas as pd


class IntegratedUiSmokeTests(unittest.TestCase):
    def _import_integrated_module(self):
        try:
            return importlib.import_module("chembl_bioactivity_integrated")
        except ModuleNotFoundError as exc:
            if exc.name not in {"chembl_webresource_client", "itables"}:
                self.skipTest(f"UI dependency unavailable: {exc.name}")

        sys.modules.pop("chembl_bioactivity_integrated", None)
        fake_chembl = types.ModuleType("chembl_webresource_client")
        fake_chembl_new_client = types.ModuleType(
            "chembl_webresource_client.new_client"
        )
        fake_chembl_new_client.new_client = types.SimpleNamespace()
        fake_itables = types.ModuleType("itables")
        fake_itables.init_notebook_mode = lambda *args, **kwargs: None
        fake_itables.show = lambda *args, **kwargs: None
        fake_itables.options = types.SimpleNamespace()
        fake_modules = {
            "chembl_webresource_client": fake_chembl,
            "chembl_webresource_client.new_client": fake_chembl_new_client,
            "itables": fake_itables,
        }
        with mock.patch.dict(sys.modules, fake_modules):
            return importlib.import_module("chembl_bioactivity_integrated")

    def test_integrated_module_imports_when_ui_dependencies_exist(self):
        module = self._import_integrated_module()

        self.assertTrue(callable(module.interactive_mode))

    def test_pubchem_basic_properties_keep_source_column_in_pubchempy_path(self):
        module = self._import_integrated_module()

        fake_pcp = mock.Mock()
        fake_pcp.get_properties.return_value = pd.DataFrame(
            [{"IUPACName": "example", "MolecularWeight": 13.0}], index=[123]
        )
        had_pcp = hasattr(module, "pcp")
        previous_pcp = getattr(module, "pcp", None)
        previous_flag = module._PCP
        try:
            module._PCP = True
            module.pcp = fake_pcp
            df = module.pubchem_basic_props_df(123)
        finally:
            module._PCP = previous_flag
            if had_pcp:
                module.pcp = previous_pcp
            else:
                delattr(module, "pcp")

        self.assertEqual(list(df.columns), ["Source", "Property", "Value"])
        self.assertEqual(set(df["Source"]), {"PubChem (computed)"})
        self.assertIn("IUPACName", set(df["Property"]))

    def test_build_activity_df_keeps_prefetched_target_name(self):
        module = self._import_integrated_module()
        activities = [
            {
                "target_chembl_id": "CHEMBL3253",
                "target_pref_name": "Albumin",
                "standard_type": "Ki",
                "standard_relation": "=",
                "standard_value": "10",
                "standard_units": "nM",
            }
        ]

        with mock.patch.object(
            module, "fetch_target_names", side_effect=AssertionError("not called")
        ):
            df = module.build_activity_df(activities)

        self.assertEqual(df.loc[0, "Target"], "Albumin")


if __name__ == "__main__":
    unittest.main()
