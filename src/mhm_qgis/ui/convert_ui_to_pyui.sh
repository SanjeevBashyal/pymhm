#!/usr/bin/env bash
# Regenerate the Qt modules in pyui/ from the Designer forms in forms/.
set -euo pipefail
cd "$(dirname "$0")"

for form in forms/*.ui; do
    pyuic5 -x "$form" -o "pyui/ui_$(basename "${form%.ui}").py"
done

pyrcc5 forms/resources.qrc -o ../resources_rc.py
