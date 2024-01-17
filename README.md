# voidDock
Docking Protocol:
    1. Use FPocket to identify binding pocket
    2. Use Scrwl4 to mutate all residues in binding pocket to Alanine
    3. Use MGLTools to convert receptor and ligand to pdbqt files
    4. Run vina docking using identified binding pocket to place box
    5  Process results into one PDB file per docking pose

# voidDock Installation

## Step 1: Create Python Environment

Create a Python environment using Conda:
```bash
## Create env
conda create -n voidDock37 python=3.7
## Activate env
conda activate voidDock37
```

## Step 2: Install MGLTools

Download MGLTools from [here](https://ccsb.scripps.edu/mgltools/downloads/).

```bash
cd /home/{username}/bin
tar -xvzf mgltools_Linux-x86_64_1.5.7.tar.gz
cd mgltools_Linux-x86_64_1.5.7
chmod +x install.sh
./install.sh
```
## Step 3: Deal with python2 requirements for MGLTools:
Add python2.7 to your PATH in your .bashrc (or .bash_proifle)
```bash
## add to .bashrc
export PATH="$PATH:/usr/bin/python2.7"
## source
source ~/.bashrc
```
Install Numpy for Python 2:

```bash
curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o get-pip.py
python2.7 get-pip.py
python2.7 -m pip install numpy
```
## Step 4: Install Python3 Libraries

Install required Python3 libraries using pip and conda:

```bash
pip3 install pytest-shutil
pip3 install pandas
conda config --add channels conda-forge
conda install fpocket
pip3 install argpass
pip3 install tqdm
pip3 install vina
```

## Step 5: Create Config File

Create a YAML script for the configuration file, e.g., `config.yaml`. Fill in the required paths:
```
dockingTargetsInfo: 
  protDir: "/home/{username}/voidDock/receptors"
  ligandDir: "/home/{username}/voidDock/ligands"
  outDir: "/home/{username}/voidDock/outputs"
  ligandOrdersCsv: "/home/{username}/voidDock/docking_commands.csv"
toolInfo:
  mglToolsDir: "/home/{username}/bin/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs"
  util24Dir: "/home/{username}/bin/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24"
```

**Note:** Replace `{username}` and update file paths accordingly in the above commands. Make sure to adjust permissions and paths based on your system configuration. Do not change the structure of the file, aka 'dockingTargetsInfo' and 'toolINfo' is needed.

## Step 6: Run voidDock
Run the voidDock script with the provided configuration file:

```bash
python voidDock.py --config config.yaml
```

happy voidDocking!
