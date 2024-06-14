
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
from typing import Tuple, Union
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
def process_vina_results_flexible_fix(outDir: Union[PathLike, str],
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
def process_vina_results(
        dockingOrder: dict,
        outDir: str,
        dockedPdbqt: str,
        receptorPdb: str) -> None:
    '''
    After a docking simulation has been run
    Use OpenBabel to split output pdbqt file into multiple pose pdb files
    The use dataframe methods to re-merge these with receptor pdb files
    '''

    protName: str = dockingOrder["protein"]
    ligandNames: list = dockingOrder["ligands"]
    # use OBabel to split up dockedPdbt into multiple pose pdbqts
    obabelCommand: list = [
        "obabel", dockedPdbqt,
        "-O", "pose.pdb",
        "-m"]
    call(obabelCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # make a dir to put all of the docked complexes
    finalPdbDir: str = p.join(outDir, "final_docked_pdbs")
    os.makedirs(finalPdbDir, exist_ok=True)

    ligandsTag: str = "_".join(ligandNames)
    protDf: pd.DataFrame = pdbUtils.pdb2df(receptorPdb)
    for file in os.listdir(outDir):
        if not p.splitext(file)[1] == ".pdb":
            continue
        if not file.startswith("pose"):
            continue
        # extract pose number from file name
        poseNumber: str = "".join(
            [char for char in p.splitext(file)[0] if char.isdigit()])
        # load posePdb and protPdb as dataframes, concat and write back to pdb
        # file
        posePdb: str = p.join(outDir, file)
        poseDf: pd.DataFrame = pdbUtils.pdb2df(posePdb)
        complexDf: pd.DataFrame = pd.concat([protDf, poseDf])
        complexPdb: str  = p.join(finalPdbDir,f"{protName}_{ligandsTag}_{poseNumber}.pdb")
        pdbUtils.df2pdb(complexDf, complexPdb)
##########################################################################


def run_vina(outDir, configFile, ligPdbqts):
    logFile = p.join(outDir, "vina_docking.log")
    ligands = " ".join(ligPdbqts)
    with open(logFile, "a") as logFile:
        run(f"vina --config {configFile} --ligand {ligands}",
            shell=True, stdout=logFile)
##########################################################################
# writes a config file for a Vina docking simulation


def write_vina_config(
        outDir,
        receptorPdbqt,
        boxCenter,
        boxSize,
        flexPdbqt=None,
        exhaustiveness=16,
        numModes=10,
        cpus=2,
        energyRange=5,
        seed=42,
        flex=False):
    vinaConfigFile = p.join(outDir, f"vina_conf.txt")

    with open(vinaConfigFile, "w") as outFile:
        if not flex:
            outFile.write(f"receptor = {receptorPdbqt}\n")
        else:
            outFile.write(f"receptor = {receptorPdbqt}\n")
            outFile.write(f"flex = {flexPdbqt}\n\n")

        outFile.write(f"center_x = {str(boxCenter[0])}\n")
        outFile.write(f"center_y = {str(boxCenter[1])}\n")
        outFile.write(f"center_z = {str(boxCenter[2])}\n\n")

        outFile.write(f"size_x = {str(boxSize)}\n")
        outFile.write(f"size_y = {str(boxSize)}\n")
        outFile.write(f"size_z = {str(boxSize)}\n\n")

        outFile.write(f"exhaustiveness = {str(exhaustiveness)}\n")
        outFile.write(f"num_modes = {str(numModes)}\n")
        outFile.write(f"energy_range = {str(energyRange)}\n\n")
        outFile.write(f"seed = {str(seed)}\n\n")

        if not flex:
            dockedPdbqt = p.join(outDir, f"binding_poses.pdbqt")
            outFile.write(f"out = {dockedPdbqt}\n")
        else:
            dockedPdbqt = p.join(outDir, f"binding_poses.pdbqt")
            outFile.write(f"out = {dockedPdbqt}\n")
        outFile.write(f"cpu = {cpus}")

        return vinaConfigFile, dockedPdbqt

def pdbDf2seq(df: pd.DataFrame) -> dict:
    ## init dict
    aminoAcidsThreeToOne = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}
    ##
    df['oneLetter'] = df['RES_NAME'].map(aminoAcidsThreeToOne)
    uniqueResidues = df[['CHAIN_ID', 'RES_ID', 'oneLetter']].drop_duplicates(subset=['CHAIN_ID', 'RES_ID'])

    chainGroups = uniqueResidues.sort_values(by=['CHAIN_ID', 'RES_ID']).groupby('CHAIN_ID')
    sequences = {}
    for chain, group in chainGroups:
        sequences[chain] = ''.join(group['oneLetter'])
    # Print all sequences, each represented by a different chain

    return sequences
##########################################################################
def pocket_residues_to_alainine(
        protName,
        pdbFile,
        residuesToAlanine,
        dockingOrder,
        outDir):
    keepResidues = dockingOrder["keepResidues"]

    # remove residues in keepResidues, these will not be changed to Alanine
    residuesToAlanine = [residue for residue in residuesToAlanine if residue not in keepResidues]
    protDf = pdbUtils.pdb2df(pdbFile)
    pocketResidues_set = set((d['CHAIN_ID'], d['RES_ID'])
                             for d in residuesToAlanine)
    ## change residue names to new sequence
    def change_res_name(row):
        if (row['CHAIN_ID'], row['RES_ID']) in pocketResidues_set:
            return 'ALA'
        else:
            return row['RES_NAME']
    protDf['RES_NAME'] = protDf.apply(change_res_name, axis=1)


    ## get dict of one-letter sequences 
    newSequences  = pdbDf2seq(protDf)


    tmpChainPdbs = {}
    for chainId in protDf['CHAIN_ID'].unique():
        chainDf = protDf[protDf['CHAIN_ID'] == chainId]
        tmpChainPdb = p.join(outDir, f"{protName}_{chainId}.pdb")
        pdbUtils.df2pdb(chainDf, tmpChainPdb)

        chainAmpal = ampal.load_pdb(tmpChainPdb)
        newChainSequence = newSequences[chainId]
        packedChain = modelling.pack_side_chains_scwrl(chainAmpal, [newChainSequence])
        finalPdbString = packedChain.make_pdb(ligands=False)
        with open(tmpChainPdb, "w") as file:
            file.write(finalPdbString)
        tmpChainPdbs[chainId] = tmpChainPdb

    tmpChainPdbList = [tmpChainPdbs[chainId] for chainId in tmpChainPdbs]
    alaPdb = p.join(outDir, f"{protName}_pocketAla.pdb")
    pdbUtils.mergePdbs(tmpChainPdbList, alaPdb)
    for pdbFile in tmpChainPdbList:
        os.remove(pdbFile)

    return alaPdb
##########################################################################


def pdb_to_pdbqt(inPdb, outDir, jobType):
    name = p.splitext(p.basename(inPdb))[0]
    outPdbqt = p.join(outDir, f"{name}.pdbqt")
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
    call(obabelCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    with open(outPdbqt, "r") as f:
        fileContents = f.read()
    fileContents = fileContents.replace("NA", "N ")

    with open(outPdbqt, "w") as f:
        f.write(fileContents)

    return outPdbqt

##########################################################################


def set_up_directory(outDir, pathInfo, dockingOrder):
    # read protein pdb file, get name and make new dir for docking, copy over
    # protein pdb
    protName = dockingOrder["protein"]
    ligands = dockingOrder["ligands"]
    ligandDir = pathInfo["ligandDir"]
    protDir = pathInfo["protDir"]

    pocketTag = ""
    if "pocketTag" in dockingOrder:
        pocketTag = dockingOrder["pocketTag"]

    protPdb = p.join(protDir, f"{protName}.pdb")

    # read ligand pdb and copy to new run directory
    ligandNames = []
    ligPdbqts = []
    for ligandName in ligands:
        ligandNames.append(ligandName)
        ligPdbqt = p.join(ligandDir, f"{ligandName}.pdbqt")
        ligPdbqts.append(ligPdbqt)

    ligandTag = "_".join(ligandNames)
    runDir = p.join(outDir, f"{protName}_{ligandTag}_{pocketTag}")
    os.makedirs(runDir, exist_ok=True)
    # copy over protPdb, move location of var
    protPdbDest = p.join(runDir, f"{protName}.pdb")
    copy(protPdb, protPdbDest)
    protPdb = protPdbDest

    return protName, protPdb, ligPdbqts, runDir
##########################################################################


def directed_fpocket(protName, runDir, pdbFile, targetPocketResidues):
    # print("----->\tRunning Fpocket!")
    os.chdir(runDir)
    minSphereSize = "3.0"
    maxSphereSize = "6.0"
    call(["fpocket", "-f", pdbFile, "-m", minSphereSize,
         "-M", maxSphereSize], stdout=subprocess.PIPE)
    fpocketOutDir = p.join(runDir, f"{protName}_out", "pockets")
    targetResidueCounts = {}
    # look at all the pockets that FPocket finds, count how many target
    # residues are in each pocket
    for file in os.listdir(fpocketOutDir):
        targetResidueCount = 0
        fileData = p.splitext(file)
        if not fileData[1] == ".pdb":
            continue
        pocketPdb = p.join(fpocketOutDir, file)
        pocketDf = pdbUtils.pdb2df(pocketPdb)
        for targetRes in targetPocketResidues:

            targetDf = pocketDf[(pocketDf["CHAIN_ID"] == targetRes["CHAIN_ID"])
                                & (pocketDf["RES_NAME"] == targetRes["RES_NAME"])
                                & (pocketDf["RES_ID"] == targetRes["RES_ID"])]
            if len(targetDf) > 0:
                targetResidueCount += 1
        pocketId = fileData[0].split("_")[0]
        targetResidueCounts.update({pocketId: targetResidueCount})

    # select the pocket with the most active site residues
    bindingPocketId = max(targetResidueCounts, key=targetResidueCounts.get)
    bindingPocketPdb = p.join(fpocketOutDir, f"{bindingPocketId}_atm.pdb")
    bindingPocketDf = pdbUtils.pdb2df(bindingPocketPdb)

    boxCenter = [
        bindingPocketDf["X"].mean(),
        bindingPocketDf["Y"].mean(),
        bindingPocketDf["Z"].mean()]

    pocketChains = bindingPocketDf["CHAIN_ID"].to_list()
    pocketResNames = bindingPocketDf["RES_NAME"].to_list()
    pocketResIds = bindingPocketDf["RES_ID"].to_list()

    uniqueResidues = set(zip(pocketChains, pocketResNames, pocketResIds))
    pocketResidues = [{"CHAIN_ID": chainId,
                       "RES_NAME": resName,
                       "RES_ID": resId} for chainId,
                      resName,
                      resId in uniqueResidues]
    pocketResidues = sorted(
        pocketResidues, key=lambda x: (
            x['CHAIN_ID'], x['RES_ID']))

    # dump to yaml
    pocketResiduesYaml = p.join(runDir, f"{protName}_pocket_residues.yaml")
    with open(pocketResiduesYaml, "w") as yamlFile:
        yaml.dump(pocketResidues, yamlFile, default_flow_style=False)

    # print("----->\tFpocket Success!")

    return boxCenter, pocketResidues
##########################################################################
