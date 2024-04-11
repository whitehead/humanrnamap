## Pre-requisites

- Python 3.11 or higher
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
   ```

## Installation

### Virtual Environment & Dependencies

1. Create a virtual environment:
    ```bash
    python -m venv venv
    ```

2. Activate the virtual environment:
    - Windows:
        ```bash
        venv\Scripts\activate
        ```
    - Linux/macOS:
        ```bash
        source venv/bin/activate
        ```

3. Install the dependencies:
    ```bash
    pip install -r requirements.txt
    ```
