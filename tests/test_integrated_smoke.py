import importlib
import unittest


class IntegratedUiSmokeTests(unittest.TestCase):
    def test_integrated_module_imports_when_ui_dependencies_exist(self):
        try:
            module = importlib.import_module("chembl_bioactivity_integrated")
        except ModuleNotFoundError as exc:
            self.skipTest(f"UI dependency unavailable: {exc.name}")

        self.assertTrue(callable(module.interactive_mode))


if __name__ == "__main__":
    unittest.main()
