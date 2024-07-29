#!/bin/bash
# Pass it a list of svg files.
# e.g. ./append_svg_images.sh main.svg legend.svg -output output.svg
set -e
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source "$SCRIPT_DIR/venv/bin/activate"
python "$SCRIPT_DIR/append_svgs.py" "$@"
