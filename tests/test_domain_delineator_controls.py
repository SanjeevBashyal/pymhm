from pathlib import Path


def test_domain_outlet_control_is_never_disabled_by_its_value():
    """An unchanged outlet location is still a valid domain/gauge location.

    The enabling lives in the form's controller; the wiring that must NOT
    re-disable it lives in the bindings module.
    """
    root = Path(__file__).parents[1] / "src" / "mhm_qgis"
    ui = (root / "qt" / "ui" / "forms" / "domain_delineator_dialog.ui").read_text(
        encoding="utf-8"
    )
    controller = (root / "qt" / "controllers" / "domain_delineator.py").read_text(
        encoding="utf-8"
    )
    bindings = (root / "qt" / "bindings" / "domain_delineator.py").read_text(
        encoding="utf-8"
    )

    assert 'name="widget_domainControls"' in ui
    assert "dialog.widget_domainControls.setEnabled(True)" in controller
    assert "checkBox_isDomainOutlet.toggled.connect" not in bindings
