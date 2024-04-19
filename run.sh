#!/bin/bash
source venv/bin/activate

export DATAPATH=$(realpath ./bin/data_tables)
export PATH="./bin:$PATH"
python ./main.py $@
