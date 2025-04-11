## Pre-requisites

- Clone this repository
    ```bash
    git clone https://github.com/whitehead/humanrnamap
    cd humanrnamap
    ```
- Download our dataset
    ```bash
    wget -O humanrnamap_dataset.tgz https://zenodo.org/records/15198110
    tar zxvf humanrnamap_dataset.tgz
    ```

- Python 3.11 or higher
- VARNA
   https://varna.lisn.upsaclay.fr/
   ```bash
   wget https://varna.lisn.upsaclay.fr/bin/VARNAv3-93.jar
   mv VARNAv3-93.jar bin/.
   ```

- RNAstructure - efn2, Fold, ct2dot
   https://rna.urmc.rochester.edu/RNAstructureDownload.html
   efn2: Version 6.4 (December 8, 2021).
   Mathews Lab, University of Rochester.

   ```bash
   wget https://rna.urmc.rochester.edu/Releases/current/RNAstructureLinuxTextInterfaces64bit.tgz
   tar xvfz RNAstructureLinuxTextInterfaces64bit.tgz
   cd RNAstructure
   make Fold
   make efn2
   make ct2dot
   cp exe/Fold ../bin/.
   cp exe/efn2 ../bin/.
   cp exe/ct2dot ../bin/.
   cp -r data_tables ../bin/.
   ```

## Installation

### Virtual Environment & Dependencies

1. Create a virtual environment:
    ```bash
    python -m venv venv
    ```

2. Activate the virtual environment:
    ```bash
    source venv/bin/activate
    ```

3. Install the dependencies:
    ```bash
    pip install -r requirements.txt
    ```
