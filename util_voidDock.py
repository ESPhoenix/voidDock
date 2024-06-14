
# import basic libraries
import os
import subprocess
from subprocess import call, run
import os.path as p
from shutil import copy, move, rmtree
import pandas as pd
import ampal
import isambard.modelling as modelling
import yaml

from pdbUtils import pdbUtils
##########################################################################


def gen_ligand_pdbqts(dockingOrders, ligandDir):
    allLigands = []
    for dockingOrder in dockingOrders:

        ligands = dockingOrder["ligands"]
        for ligand in ligands:
            allLigands.append(ligand)

    allLigands = list(set(allLigands))
    for ligand in allLigands:
        ligPdb = p.join(ligandDir, f"{ligand}.pdb")
        if not p.isfile(ligPdb):
            print(f"{ligPdb} not found, skipping...")
            continue
        pdb_to_pdbqt(ligPdb, ligandDir, jobType="ligand")
##########################################################################


def collate_docked_pdbs(outDir, rmDirs=True):
    for dir in os.listdir(outDir):

        runDir = p.join(outDir, dir)
        finalDockedPdbDir = p.join(runDir, "final_docked_pdbs")
        if not p.isdir(finalDockedPdbDir):
            continue
        for file in os.listdir(finalDockedPdbDir):
            if not p.splitext(file)[1] == ".pdb":
                continue
            pdbFile = p.join(finalDockedPdbDir, file)
            pdbDest = p.join(outDir, file)
            copy(pdbFile, pdbDest)
        for file in os.listdir(runDir):
            if file.endswith("_pocket_residues.yaml"):
                pocketResiduesYaml = p.join(runDir, file)
                pocketResiduesDest = p.join(outDir, file)
                copy(pocketResiduesYaml, pocketResiduesDest)

        if rmDirs:
            rmtree(runDir)
##########################################################################
def gen_flex_pdbqts(protPdb, flexibleResidues, outDir):
    name = p.splitext(p.basename(protPdb))[0]
    # load pdbfile into df
    protDf = pdbUtils.pdb2df(protPdb)
    flexIndexes = []
    dfsToConcat = []
    for residue in flexibleResidues:

        chainDf = protDf[(protDf["CHAIN_ID"] == residue["CHAIN_ID"]) & 
                         (protDf["RES_ID"] == residue["RES_ID"]) &
                         (~protDf["ATOM_NAME"].isin(["CA","C","O","N"]))]

        chainIndexes = chainDf.index.to_list()
        flexIndexes += chainIndexes
        dfsToConcat.append(chainDf)
    
    flexDf = pd.concat(dfsToConcat, axis = 0)
    rigidDf = protDf.drop(index=flexIndexes)

    flexPdb = p.join(outDir,f"{name}_flex.pdb")
    rigidPdb = p.join(outDir,f"{name}_rigid.pdb")
    pdbUtils.df2pdb(flexDf,flexPdb)
    pdbUtils.df2pdb(rigidDf,rigidPdb)
    pdb_to_pdbqt(flexPdb, outDir, jobType="flex")
    pdb_to_pdbqt(rigidPdb, outDir, jobType="rigid")
    flexPdbqt = p.join(outDir,f"{name}_flex.pdbqt")
    rigidPdbqt = p.join(outDir,f"{name}_rigid.pdbqt")

    return rigidPdbqt, flexPdbqt


##########################################################################
def process_vina_results(
        dockingOrder,
        outDir,
        dockedPdbqt,
        receptorPdb):

    protName = dockingOrder["protein"]
    ligandNames = dockingOrder["ligands"]
    # use OBabel to split up dockedPdbt into multiple pose pdbqts
    obabelCommand = [
        "obabel", dockedPdbqt,
        "-O", "pose.pdb",
        "-m"]
    call(obabelCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # make a dir to put all of the docked complexes
    finalPdbDir = p.join(outDir, "final_docked_pdbs")
    os.makedirs(finalPdbDir, exist_ok=True)

    ligandsTag = "_".join(ligandNames)
    protDf = pdbUtils.pdb2df(receptorPdb)
    for file in os.listdir(outDir):
        if not p.splitext(file)[1] == ".pdb":
            continue
        if not file.startswith("pose"):
            continue
        # extract pose number from file name
        poseNumber = "".join(
            [char for char in p.splitext(file)[0] if char.isdigit()])
        # load posePdb and protPdb as dataframes, concat and write back to pdb
        # file
        posePdb = p.join(outDir, file)
        poseDf = pdbUtils.pdb2df(posePdb)
        complexDf = pd.concat([protDf, poseDf])
        complexPdb = p.join(
            finalPdbDir,
            f"{protName}_{ligandsTag}_{poseNumber}.pdb")
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
