# voidDock
Docking Protocol:
    1. Use FPocket to identify binding pocket
    2. Use Scrwl4 to mutate all residues in binding pocket to Alanine
    3. Use MGLTools to convert receptor and ligand to pdbqt files
    4. Run vina docking using identified binding pocket to place box
    5  Process results into one PDB file per docking pose

# voidDock Installation
## Step 1: Clone this repository
```bash
git clone https://github.com/ESPhoenix/voidDock
```

## Step 2: Create Python Environment

Create a Python environment using Conda:
```bash
## Create env
conda create -n voidDock37 python=3.7
## Activate env
conda activate voidDock37
```

## Step 3: Install Python3 Libraries

Install required Python libraries:

```bash
pip install -r requirements.txt
```

## Step 4: Create Config File

Create a YAML script for the configuration file, e.g., `config.yaml`. Fill in the required paths and information:
```yaml
#### protDir    : path to directory containing proten (receptor) pdb files
#### ligandDir  : path to directory containing ligand pdb files  
#### outputDir  : path of desired output directory (does not need to exist yet)
#### 
dockingTargetsInfo: 
  protDir: "/home/{username}/voidDock/receptors"
  ligandDir: "/home/{username}/voidDock/ligands"
  outDir: "/home/{username}/voidDock/outputs"

#### totalCpuUsage  :   number of cpu cores the whole script will use
#### cpusPerRun     :   number of cpu cores each vina process will use
#### total vina processes at a time will be:
#### totalCpuUsage // cpusPerRun
####
cpuInfo:
  totalCpuUsage: 8
  cpusPerRun: 2

#### list of dictionaries containing docking information for individual docking runs
#### each dictionary contains:
#### protein  :   protein name (this is the pdb file name)
#### ligand   :   ligand name (this is the pdb file name)
#### pocketResidues :   list of residues that will help identify the binding pocket in the format "CHAIN_ID:RES_NAME:RES_ID"
####
dockingOrders:
  - protein: "134189607" 
    ligand: "Lumoflavin"
    pocketResidues: ["A:GLY:371", "A:ALA:854", "A:ALA:584", "A:GLY:171"]
  - protein: "134189607"  
    ligand: "FAD"
    pocketResidues: ["A:GLY:371", "A:ALA:854", "A:ALA:584", "A:GLY:171"]
  - protein: "62562582"
    ligand: "Lumoflavin"
    pocketResidues: ["A:ASN:60", "A:LYS:78", "A:ALA:77", "A:ILE:62"]
```

**Note:** Replace `{username}` and update file paths accordingly in the above commands.
Replace all other information with your desired inputs.
Make sure to adjust permissions and paths based on your system configuration. Do not change the structure of the file, aka 'dockingTargetsInfo' 

## Step 5: Run voidDock
Run the voidDock script with the provided configuration file:

```bash
python voidDock.py --config config.yaml
```

happy voidDocking!
