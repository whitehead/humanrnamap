#!/bin/bash
# Pass it a list of svg files.
# e.g. ./trim_svg_images.sh *.svg
set -e
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source "$SCRIPT_DIR/venv/bin/activate"
python "$SCRIPT_DIR/trim_svgs.py" "$@"
