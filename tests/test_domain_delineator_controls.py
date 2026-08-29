from pathlib import Path


def test_domain_outlet_control_is_never_disabled_by_its_value():
    root = Path(__file__).parents[1] / "src" / "mhm_qgis"
    ui = (root / "ui" / "domain_delineator_dialog.ui").read_text(encoding="utf-8")
    controller = (root / "domain_delineator_dialog.py").read_text(encoding="utf-8")

    assert 'name="widget_domainControls"' in ui
    assert "self.widget_domainControls.setEnabled(True)" in controller
    assert "self.checkBox_isDomainOutlet.toggled.connect(\n            self.widget_domainControls.setEnabled)" not in controller
