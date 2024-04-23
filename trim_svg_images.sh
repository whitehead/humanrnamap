#!/bin/bash
# Pass it a list of svg files.
# e.g. ./trim_svg_images.sh *.svg
set -e
source venv/bin/activate
python ./trim_svgs.py $@
