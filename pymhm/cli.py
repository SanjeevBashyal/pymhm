#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Command line interface for PymHM."""

import argparse
import sys


def launch_gui():
    """Launch the plugin dialog with standalone file-path widgets."""
    from .standalone_qgis import install

    install(force=True)

    from qgis.PyQt.QtWidgets import QApplication

    app = QApplication.instance()
    owns_app = app is None
    if app is None:
        app = QApplication(sys.argv[:1])

    from .pymhm_dialog import pymhmDialog

    dialog = pymhmDialog()
    dialog.setWindowTitle("PymHM")
    dialog.show()
    if owns_app:
        return app.exec_()
    return dialog.exec_()


def main(argv=None):
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="PymHM - Python Mesoscale Hydrological Model",
        prog="pymhm"
    )
    
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s 0.1.0"
    )
    
    parser.add_argument(
        "--info",
        action="store_true",
        help="Show package information"
    )

    parser.add_argument(
        "--no-gui",
        action="store_true",
        help="Do not open the standalone graphical interface."
    )

    args = parser.parse_args(argv)
    
    if args.info:
        print("PymHM - Python Mesoscale Hydrological Model")
        print("Version: 0.1.0")
        print("Author: Sanjeev Bashyal")
        print("Email: sanjeev.bashyal01@gmail.com")
        print("Description: Python package for mesoscale Hydrological Model")
        print("Homepage: https://github.com/SanjeevBashyal/pymhm")
        return 0

    if args.no_gui:
        print("PymHM CLI installed. Run `pymhm` without --no-gui to open the GUI.")
        return 0
    
    try:
        return launch_gui()
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
