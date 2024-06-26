
# import basic libraries
import os
import subprocess
from subprocess import call, run
import os.path as p
from shutil import copy, move, rmtree
import pandas as pd
import yaml
## WellsWood Lab libraries
import ampal
import isambard.modelling as modelling
from pdbUtils import pdbUtils
## clean code
from typing import Tuple, Union, Literal
from os import PathLike
import numpy as np
##########################################################################

def gen_ligand_pdbqts(dockingOrders: dict, ligandDir: Union[PathLike, str]) -> None:
    '''
    Before running docking, generate pdbqt files for all ligands
    This saves us from generating a ligand pdbqt per docking run
    '''
    ## loop through ligand pdb files in the docking orders
    ## append these files to a list
    ## get unique entries
    allLigands: list = []
    for dockingOrder in dockingOrders:
        ligands: list = dockingOrder["ligands"]
        for ligand in ligands:
            allLigands.append(ligand)
    allLigands: list = list(set(allLigands))
    ## look in ligandsDir for these ligands
    ## call pdb_to_pdbqt to convert to pdbqt file
    for ligand in allLigands:
        ligPdb: Union[PathLike, str] = p.join(ligandDir, f"{ligand}.pdb")
        if not p.isfile(ligPdb):
            print(f"{ligPdb} not found, skipping...")
            continue
        pdb_to_pdbqt(ligPdb, ligandDir, jobType="ligand")
##########################################################################


def collate_docked_pdbs(outDir: Union[PathLike, str], rmDirs: bool = True) -> None:
    '''
    After a docking simulation has been run, get:
     - output pdb files
     - pocket_residues.yaml (useful for further pipelines)
     Move these files to a new dir
     Remove the rest of the outputs to save space
    '''
    ## loop through run directory and copy required files to final_docked_pdbs
    for dir in os.listdir(outDir):
        runDir: Union[PathLike, str]= p.join(outDir, dir)
        finalDockedPdbDir: Union[PathLike, str]= p.join(runDir, "final_docked_pdbs")
        if not p.isdir(finalDockedPdbDir):
            continue
        for file in os.listdir(finalDockedPdbDir):
            if not p.splitext(file)[1] == ".pdb":
                continue
            pdbFile: Union[PathLike, str]= p.join(finalDockedPdbDir, file)
            pdbDest: Union[PathLike, str]= p.join(outDir, file)
            copy(pdbFile, pdbDest)
        for file in os.listdir(runDir):
            if file.endswith("_pocket_residues.yaml"):
                pocketResiduesYaml: Union[PathLike, str]= p.join(runDir, file)
                pocketResiduesDest: Union[PathLike, str]= p.join(outDir, file)
                copy(pocketResiduesYaml, pocketResiduesDest)
        ## if specified in config, remove all other files and directories
        ## produced by fpocket/docking procedure
        if rmDirs:
            rmtree(runDir)

##########################################################################
def get_ligand_res_names(ligandDir: os.PathLike, ligands: list) -> list:
    """
    Simple loop function to look through ligand pdb files and extract three-letter residue names
    """
    resNames: list = []
    for ligand in ligands:
        ligandPdb: Union[PathLike, str] = p.join(ligandDir, f"{ligand}.pdb")
        ligandDf: pd.DataFrame = pdbUtils.pdb2df(ligandPdb)
        resName: str = ligandDf["RES_NAME"].tolist()[0]
        resNames.append(resName)

    return resNames
    





##########################################################################
def read_docking_results(dockedPdbqt: Union[PathLike, str]) -> list:
    """
    Reads docking result pdbqt file into one dataframe per binding pose
    Returns these dataframes in a list
    """
    # remove ROOT/BRANCH
    pdbqtColumns: list    =   ["ATOM","ATOM_ID", "ATOM_NAME", "RES_NAME",
                    "CHAIN_ID", "RES_ID", "X", "Y", "Z", "OCCUPANCY", 
                    "BETAFACTOR","CHARGE"]
    columsNums: list[tuple] = [(0, 6), (6, 11), (11, 17), (17, 21), (21, 22), (22, 26), 
                  (26, 38), (38, 46), (46, 54), (54, 60), (60, 70), (70, 77)]
    # read pdbqt file into multiple dataframes
    dockingDfList: list =[]
    data: list = []
    # read output PDBQT file into set of dataframes
    with open(dockedPdbqt, 'r') as file:            
        for line in file:   
            if line.startswith("MODEL"):        # Each binding pose starts with "MODEL"
                if data == []:                  # skip 1st "MODEL"
                    continue
                df: pd.DataFrame = pd.DataFrame(data,columns=pdbqtColumns)
                df[["ATOM_ID","RES_ID"]] = df[["ATOM_ID","RES_ID"]].astype(int)
                df[["X","Y","Z","OCCUPANCY","BETAFACTOR"]]=df[["X","Y","Z","OCCUPANCY","BETAFACTOR"]].astype(float)
                dockingDfList.append(df)
                data = []
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                record: list = [line[start:end].strip() for start, end in columsNums]
                data.append(record)
    # deal with last entry in pdbqtfile
    df: pd.DataFrame = pd.DataFrame(data,columns=pdbqtColumns)
    df[["ATOM_ID","RES_ID"]] = df[["ATOM_ID","RES_ID"]].astype(int)
    df[["X","Y","Z","OCCUPANCY","BETAFACTOR"]]=df[["X","Y","Z","OCCUPANCY","BETAFACTOR"]].astype(float)
    dockingDfList.append(df)

    return dockingDfList


##########################################################################

def gen_ligand_pdbqts(dockingOrders: dict, ligandDir: Union[PathLike, str]) -> None:
    '''
    Before running docking, generate pdbqt files for all ligands
    This saves us from generating a ligand pdbqt per docking run
    '''
    ## loop through ligand pdb files in the docking orders
    ## append these files to a list
    ## get unique entries
    allLigands: list = []
    for dockingOrder in dockingOrders:
        ligands: list = dockingOrder["ligands"]
        for ligand in ligands:
            allLigands.append(ligand)
    allLigands: list = list(set(allLigands))
    ## look in ligandsDir for these ligands
    ## call pdb_to_pdbqt to convert to pdbqt file
    for ligand in allLigands:
        ligPdb: Union[PathLike, str] = p.join(ligandDir, f"{ligand}.pdb")
        if not p.isfile(ligPdb):
            print(f"{ligPdb} not found, skipping...")
            continue
        pdb_to_pdbqt(ligPdb, ligandDir, jobType="ligand")
##########################################################################

def gen_flex_pdbqts(protPdb: Union[PathLike, str],
                     flexibleResidues: list,
                       outDir: Union[PathLike, str]) -> Tuple[str, str]:
    """
    Using user generated list of desired flexible residues:
    1. Create separate pdb files for flexible residues
    2. Create new rigid pdb file by dropping flexible residues from protein pdb
    3. Convert these pdb files to pdbqt files
    """

    ## get protein name from its filename
    protName: Union[PathLike, str] = p.splitext(p.basename(protPdb))[0]
    # load pdbfile into df
    protDf: pd.DataFrame = pdbUtils.pdb2df(protPdb)
    ## create empty lists to store data
    flexIndexes: list = []
    dfsToConcat: list = []
    ## use flexibleResidues entry in config to get dataframes of sidechains
    for residue in flexibleResidues:
        chainDf: pd.DataFrame = protDf[(protDf["CHAIN_ID"] == residue["CHAIN_ID"]) & 
                         (protDf["RES_ID"] == residue["RES_ID"]) &
                         (~protDf["ATOM_NAME"].isin(["CA","C","O","N"]))]

        chainIndexes: list = chainDf.index.to_list()
        flexIndexes += chainIndexes
        dfsToConcat.append(chainDf)
    ## concat flexible sidechain dataframes, drop from protDf to create rigidDf
    flexDf: pd.DataFrame = pd.concat(dfsToConcat, axis = 0)
    rigidDf: pd.DataFrame = protDf.drop(index=flexIndexes)
    ## write dataframes to pdb files
    flexPdb: Union[PathLike, str] = p.join(outDir,f"{protName}_flex.pdb")
    rigidPdb: Union[PathLike, str] = p.join(outDir,f"{protName}_rigid.pdb")
    pdbUtils.df2pdb(flexDf,flexPdb)
    pdbUtils.df2pdb(rigidDf,rigidPdb)
    ## convert pdb files to pdbqt files
    pdb_to_pdbqt(flexPdb, outDir, jobType="flex")
    pdb_to_pdbqt(rigidPdb, outDir, jobType="rigid")
    flexPdbqt: Union[PathLike, str] = p.join(outDir,f"{protName}_flex.pdbqt")
    rigidPdbqt: Union[PathLike, str] = p.join(outDir,f"{protName}_rigid.pdbqt")

    return rigidPdbqt, flexPdbqt

##########################################################################
def process_vina_results(outDir: Union[PathLike, str],
                                        dockedPdbqt: Union[PathLike, str],
                                        receptorPdbqt: Union[PathLike, str],
                                        dockingOrder: dict,
                                        ligandDir: Union[PathLike, str]) -> list:
    """
    Main function for processing vina outputs
    1. Reads binding_poses.pdbqt into one DataFrame per pose
    2. Merges individual poses with rigid receptor 
    3. Reorganises flexible residues back into correct place
    4. Ensures that ligands have a different chain Id to the receptor
    5. Writes processed poses to pdb files
    """

    # read output pdbqt file into a list of dataframes
    dockingDfList: list = read_docking_results(dockedPdbqt)
    receptorDf: pd.DataFrame = pdbUtils.pdbqt2df(receptorPdbqt)
    protName: str = dockingOrder["protein"]
    ligandNames: list[str] = dockingOrder["ligands"]
    ligandResNames: list[str] = get_ligand_res_names(ligandDir, ligandNames)
    # read ligand pdb and copy to new run directory
    ligTag: str = "_".join(ligandNames)
    nameTag: str = f"{protName}_{ligTag}"
    dockedPdbs: list = splice_docking_results(dockingDfList, receptorDf, outDir, nameTag, dockingOrder, ligandResNames)
    return dockedPdbs

##########################################################################
def splice_docking_results(dockingDfList: list,
                            receptorDf: pd.DataFrame,
                              outDir: Union[PathLike, str],
                                nameTag: str,
                                  dockingOrder: dict,
                                    ligandResNames: list) -> list:
    '''
    Manual method using dataframes to splice docking dataframes with receptor dataframes
    This method must:
        1. Work for flexible residues (OpenBabel fails to do this!)
        2. Work for multiple ligands
        3. Preserve Chain Ids for the receptor
        4. Ensure that Chain Ids for ligands are different to receptor chain Ids
    '''
    ## set elements in receptor to None (all messed up from vina outputs)
    receptorDf.loc[:,"ELEMENT"] = None

    ## make an output dir to put files in 
    finalPdbDir: Union[PathLike, str] = p.join(outDir,"final_docked_pdbs")
    os.makedirs(finalPdbDir,exist_ok=True)
    ## init an empty list to store docked pdb files
    dockedPdbs: list = []
    ## loop over each pose in dockingDfList
    for poseNumber, dockedDf in zip(range(1,len(dockingDfList)+1),dockingDfList):
        dockedDf.loc[:,"ELEMENT"] = None

        ## split dockingDf into ligand and flexible residues
        ligandDf: pd.DataFrame = dockedDf[dockedDf["RES_NAME"].isin(ligandResNames)]
        flexibleResidueDf: pd.DataFrame = dockedDf[~dockedDf["RES_NAME"].isin(ligandResNames)]
         
        ## find  max chain ID in receptor, set ligand to one more than that
        chainIdIdCounter: str = sorted(receptorDf["CHAIN_ID"].unique().tolist(), reverse=True)[0]
        for ligandResName in ligandResNames:
            thisLigandIndexes: pd.Series = ligandDf["RES_NAME"] == ligandResName
            chainIdIdCounter: str = chr((ord(chainIdIdCounter) - ord('A') + 1) % 26 + ord('A'))
            ligandDf.loc[thisLigandIndexes,"CHAIN_ID"] = chainIdIdCounter
            ligandDf.loc[thisLigandIndexes, "RES_ID"] = 1
        # Concat docked and rigid DFs togeter - this is in a weird order
        wholeDisorderedDf: pd.DataFrame = pd.concat([receptorDf, flexibleResidueDf, ligandDf],axis=0)

        dfsToConcat: str = []
        for chainId in sorted(pd.unique(wholeDisorderedDf["CHAIN_ID"]).tolist()):
            chainDf: pd.DataFrame = wholeDisorderedDf[wholeDisorderedDf["CHAIN_ID"] == chainId]
            for resId in sorted(chainDf["RES_ID"].unique().tolist()):
                resDf: pd.DataFrame = chainDf[chainDf["RES_ID"] == resId]

                dfsToConcat.append(resDf)
        # concat into correct order
        wholeDf: pd.DataFrame = pd.concat(dfsToConcat)
        # re-do atom numbers
        wholeDf.loc[:,"ATOM_ID"] = range(1,len(wholeDf)+1)
        # save as pdb file
        saveFile: Union[PathLike, str] = p.join(finalPdbDir, f"{nameTag}_{str(poseNumber)}.pdb")
        pdbUtils.df2pdb(df=wholeDf,outFile=saveFile)
        dockedPdbs.append(saveFile)
    return dockedPdbs

##########################################################################


def run_vina(outDir: Union[os.PathLike, str],
              configFile: dict,
                ligPdbqts: list) -> None:
    """
    Small function for running vina docking
    We pass a space delimited string of ligand pdbqt paths to vina
    so that multiple ligands can be docked at the same time
    this does not work with the config for some reason
    """
    logFile = p.join(outDir, "vina_docking.log")
    ligands = " ".join(ligPdbqts)
    with open(logFile, "a") as logFile:
        run(f"vina --config {configFile} --ligand {ligands}",
            shell=True, stdout=logFile) ## TODO: make this work without shell=True
        
##########################################################################

def write_vina_config(
        outDir: Union[os.PathLike, str],
        receptorPdbqt: Union[os.PathLike, str],
        boxCenter: list,
        boxSize: list,
        flexPdbqt: Union[os.PathLike, str]=None,
        exhaustiveness: int=16,
        numModes: int=10,
        cpus: int=2,
        energyRange: int=5,
        seed: int=42,
        flex: bool=False) -> Tuple[Union[os.PathLike, str],  Union[os.PathLike, str]]:
    """
    Writes a config file to be passed to vina
    We have selected some sensible defaults, but you can modify these to your liking 
    using the voidDock config.yaml file
    """
    vinaConfigFile: Union[os.PathLike, str] = p.join(outDir, f"vina_conf.txt")

    with open(vinaConfigFile, "w") as outFile:
        ## input pdbqt files, use flexible residues if required
        outFile.write(f"receptor = {receptorPdbqt}\n")
        if flex:
            outFile.write(f"flex = {flexPdbqt}\n\n")
        ## docking box center coords (calculated as center of pocket)
        outFile.write(f"center_x = {str(boxCenter[0])}\n")
        outFile.write(f"center_y = {str(boxCenter[1])}\n")
        outFile.write(f"center_z = {str(boxCenter[2])}\n\n")
        ## docking box size
        outFile.write(f"size_x = {str(boxSize)}\n")
        outFile.write(f"size_y = {str(boxSize)}\n")
        outFile.write(f"size_z = {str(boxSize)}\n\n")
        ## vina parameters
        outFile.write(f"exhaustiveness = {str(exhaustiveness)}\n")
        outFile.write(f"num_modes = {str(numModes)}\n")
        outFile.write(f"energy_range = {str(energyRange)}\n\n")
        outFile.write(f"seed = {str(seed)}\n\n")
        outFile.write(f"cpu = {cpus}\n\n")
        ## output pdbqt file
        dockedPdbqt = p.join(outDir, "binding_poses.pdbqt")
        outFile.write(f"out = {dockedPdbqt}\n")

        return vinaConfigFile, dockedPdbqt
##########################################################################

def pdb_to_pdbqt(inPdb: Union[os.PathLike, str],
                outDir: Union[os.PathLike, str],
                jobType: Literal["flex","rigid","ligand"]) -> Union[os.PathLike, str]:
    """
    Uses openBabel to convert pdb files to pdbqt files for vina inputs
    Accepts "flex", "rigid", "ligand" jobTypes depending on what degree of flexibility is required
    """
    ## get name from input pdb file, use that to name new pdbqt file
    name = p.splitext(p.basename(inPdb))[0]
    outPdbqt = p.join(outDir, f"{name}.pdbqt")

    ## depending on jobtype, create an openBabel command
    ## for flexible residues (the -s flag is special here)
    if jobType == "flex":
        obabelCommand = [
            "obabel",
            "-i",
            "pdb",
            inPdb,
            "-o",
            "pdbqt",
            "-O",
            outPdbqt,
            "-xs"]
    ## for rigid receptors (the -r flag is special here)
    elif jobType == "rigid":
        obabelCommand = [
            "obabel",
            "-i",
            "pdb",
            inPdb,
            "-o",
            "pdbqt",
            "-O",
            outPdbqt,
            "-xr"]
    ## for rigid receptors (the -r flag is special here)
    elif jobType == "ligand":
        obabelCommand = [
            "obabel",
            "-i",
            "pdb",
            inPdb,
            "-o",
            "pdbqt",
            "-O",
            outPdbqt,
            "-xn"]
    ## call the openBabel command 
    call(obabelCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ## openBabel sometimes makes Nitrogen atoms into sodiums, quick fix here!
    with open(outPdbqt, "r") as f:
        fileContents = f.read()
    fileContents = fileContents.replace("NA", "N ")

    with open(outPdbqt, "w") as f:
        f.write(fileContents)

    return outPdbqt

##########################################################################


def set_up_directory(outDir: Union[os.PathLike, str],
                    pathInfo: dict,
                    dockingOrder: dict) -> Tuple[str, Union[os.PathLike, str], list, Union[os.PathLike, str]]:
    """
    Set up a run directory for an individual docking simulation:
    1. Read pathInfo and dockingOrder dictionaries to get locations of input files
    2. Find receptor pdb file and ligand pdbqt files 
            (which have been created prior to running any docking simulations)
    3. Create a new directory to perform this docking simulation
    4. Copy protein pdb file over to this new diretory

    Returns:
    - protName (str): The name of the receptor (taken from the pdb file)
    - protPdb (PATH): The full path to the receptor pdb file
    - ligPdbqts (list): A list containing full paths to ligand pdbqt files
    - runDir (PATH): The full path to the directory made in this function
    """

    # read protein pdb file, get name and make new dir for docking, copy over
    # protein pdb
    protName: str = dockingOrder["protein"]
    ligands: list = dockingOrder["ligands"]
    ligandDir: Union[os.PathLike, str] = pathInfo["ligandDir"]
    protDir: Union[os.PathLike, str] = pathInfo["protDir"]

    protPdb: Union[os.PathLike, str] = p.join(protDir, f"{protName}.pdb")

    # read ligand pdb and copy to new run directory
    ligandNames: list = []
    ligPdbqts:  list = []
    for ligandName in ligands:
        ligandNames.append(ligandName)
        ligPdbqt: Union[os.PathLike, str] = p.join(ligandDir, f"{ligandName}.pdbqt")
        ligPdbqts.append(ligPdbqt)
    ## make a unique name for this docking simulation, 
    ## create directory with this name within outDir
    ligandTag: str = "_".join(ligandNames)
    runDir: Union[os.PathLike, str] = p.join(outDir, f"{protName}_{ligandTag}")
    os.makedirs(runDir, exist_ok=True)
    # copy over protPdb, move location of var
    protPdbDest: Union[os.PathLike, str] = p.join(runDir, f"{protName}.pdb")
    copy(protPdb, protPdbDest)
    protPdb: Union[os.PathLike, str] = protPdbDest

    return protName, protPdb, ligPdbqts, runDir

