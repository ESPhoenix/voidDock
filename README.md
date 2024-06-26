# voidDock
Docking Protocol:
    1. Use FPocket to identify binding pocket
    2. Use Scrwl4 to mutate all residues in binding pocket to Alanine
    3. Use MGLTools to convert receptor and ligand to pdbqt files
    4. Run vina docking using identified binding pocket to place box
    5  Process results into one PDB file per docking pose

# voidDock Installation
# Step 1: Clone this repository
```bash
git clone https://github.com/ESPhoenix/voidDock
```

# Step 2: Create Python Environment

Create a Python environment using Conda:
```bash
# Create env
conda create -n voidDock38 python=3.8
# Activate env
conda activate voidDock38
```

# Step 3: Install Python3 Libraries

Install required Python libraries:

```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install openbabel
conda install fpocket
conda install vina
# if this fails use: conda install conda-forge::vina
# if this fails, use sudo apt-get install autodock-vina
pip install pandas
pip install numpy
pip install pyyaml
pip install argpass
pip install cython  # must install cython before Ampal and Isambard!
pip install ampal
pip install isambard
pip install pdbutils
```
Note that if some of these installs fail, you may need to manually debug the installation(s) that failed

# Step 4: Install Scwrl4

To get the Scrwl4 program, you will need to obtain a liscense from:
[The Dunbrack Lab](http://dunbrack.fccc.edu/lab/scwrl, "link to Scwrl4 install instructions")
Once you have done this, you will be sent a link to a download location for Scrwl4
Follow the instructions in the installer
We reccomend that you do not use the default install location. Instead use ~/bin/Scrwl4 or an equivalent
After this, remember to add your install location to your PATH variable in .bashrc / .bash_profile



# Step 4: Create Config File

Create a YAML script for the configuration file, e.g., `config.yaml`. Fill in the required paths and information:
```yaml
# protDir    : path to directory containing proten (receptor) pdb files
# ligandDir  : path to directory containing ligand pdb files
# outputDir  : path of desired output directory (does not need to exist yet)
####
pathInfo:
    protDir: "/home/{username}/voidDock/receptors"
    ligandDir: "/home/{username}/voidDock/ligands"
    outDir: "/home/{username}/voidDock/outputs"

# totalCpuUsage  :   number of cpu cores the whole script will use
# cpusPerRun     :   number of cpu cores each vina process will use
# total vina processes at a time will be:
# totalCpuUsage // cpusPerRun
####
cpuInfo:
    totalCpuUsage: 8
    cpusPerRun: 2

# list of dictionaries containing docking information for individual docking runs
# each dictionary contains:
# protein           :   protein name (this is the pdb file name without extension)
# ligands           :   list of ligand names (these are the pdb file names without extensions)
# pocketResidues    :   list of residues that will help identify the binding pocket
#                       each residue is a dict containing CHAIN_ID:{str}, RES_NAME:{str}, RES_ID:{int}
# keepResidues      :   list of residues that must be conserved (this could be important for catalysis/binding)
#                       each residue is a dict containing CHAIN_ID:{str}, RES_NAME:{str}, RES_ID:{int}
####
dockingOrders:
  - protein: "A0A0D2XFD3"
    ligands: ["FFF", "TPA"]
    pocketResidues: 
      - CHAIN_ID: "A"
        RES_NAME: "PHE"
        RES_ID: 216
      - CHAIN_ID: "A"
        RES_NAME: "ARG"
        RES_ID: 217
      - CHAIN_ID: "A"
        RES_NAME: "TRP"
        RES_ID: 213
      - CHAIN_ID: "A"
        RES_NAME: "SER"
        RES_ID: 214
        
    keepResidues: 
      - CHAIN_ID: "A"
        RES_NAME: "PHE"
        RES_ID: 216
      - CHAIN_ID: "A"
        RES_NAME: "ARG"
        RES_ID: 217
      - CHAIN_ID: "A"
        RES_NAME: "TRP"
        RES_ID: 213
      - CHAIN_ID: "A"
        RES_NAME: "SER"
        RES_ID: 214
```

**Note: ** Replace `{username}` and update file paths accordingly in the above commands.
Replace all other information with your desired inputs.
Make sure to adjust permissions and paths based on your system configuration. Do not change the structure of the file, aka 'dockingTargetsInfo'

# Step 5: Run voidDock
Run the voidDock script with the provided configuration file:

```bash
python voidDock.py - -config config.yaml
```

happy voidDocking!
