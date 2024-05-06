
## import basic libraries
from itertools import chain
import os
import subprocess
from subprocess import call
import os.path as p
from shutil import copy, move
import pandas as pd
import json
import ampal 
import isambard.specifications as specifications 
import isambard.modelling as modelling
import yaml

from pdbUtils import pdbUtils

#########################################################################################################################
def process_vina_results(protName, ligandName, outDir,dockedPdbqt,receptorPdbqt):
    # read output pdbqt file into a list of dataframes
    dockingDfList = read_docking_results(dockedPdbqt)
    receptorDf = pdbUtils.pdbqt2df(receptorPdbqt)
    
    splice_docking_results(protName, ligandName, dockingDfList, receptorDf,outDir)
#########################################################################################################################
def splice_docking_results(protName, ligandName, dockingDfList, receptorDf, outDir):
    finalPdbDir = p.join(outDir,"final_docked_pdbs")
    os.makedirs(finalPdbDir,exist_ok=True)
    ## loop over each pose in dockingDfList
    for poseNumber, dockedDf in zip(range(1,len(dockingDfList)+1),dockingDfList):
        ## find  max chain ID in receptor, set ligand to one more than that
        lastChainIdInProt = receptorDf.iloc[-1]["CHAIN_ID"]
        ligandChainId = chr((ord(lastChainIdInProt) - ord('A') + 1) % 26 + ord('A'))
        dockedDf.loc[dockedDf["RES_ID"] == 0,"CHAIN_ID"] = ligandChainId
        ## Set ligand residue number to 1
        ligandResidueId = 1
        dockedDf.loc[dockedDf["RES_ID"] == 0,"RES_ID"] = ligandResidueId
        ## chage HETATM to ATOM for dockedDf
        dockedDf.loc[:,"ATOM"] = "ATOM"
        # Concat docked and rigid DFs togeter - this is in a weird order
        wholeDisorderedDf = pd.concat([dockedDf,receptorDf],axis=0)
        # get a list of unique residue Ids
        uniqueResidues = sorted(pd.unique(wholeDisorderedDf["RES_ID"]).tolist())
        # get one df per residue
        orderedResidues = []
        for residueNum in uniqueResidues:
            residueDf = wholeDisorderedDf[wholeDisorderedDf["RES_ID"]==residueNum]
            orderedResidues.append(residueDf)
        # concat into correct order
        wholeDf = pd.concat(orderedResidues)
        # re-do atom numbers
        wholeDf.loc[:,"ATOM_ID"] = range(1,len(wholeDf)+1)
        # save as pdb file
        saveFile = p.join(finalPdbDir, f"{protName}_{ligandName}_{str(poseNumber)}.pdb")
        pdbUtils.df2pdb(df=wholeDf,outFile=saveFile)

#########################################################################################################################
def read_docking_results(dockedPdbqt):
    # remove ROOT/BRANCH
    pdbqtColumns    =   ["ATOM","ATOM_ID", "ATOM_NAME", "RES_NAME",
                    "CHAIN_ID", "RES_ID", "X", "Y", "Z", "OCCUPANCY", 
                    "BETAFACTOR","CHARGE", "ELEMENT"]
    columsNums = [(0, 6), (6, 11), (11, 17), (17, 21), (21, 22), (22, 26), 
                  (26, 38), (38, 46), (46, 54), (54, 60), (60, 70), (70, 77), (77, 79)]
    # read pdbqt file into multiple dataframes
    dockingDfList =[]
    data = []
    # read output PDBQT file into set of dataframes
    with open(dockedPdbqt, 'r') as file:            
        for line in file:   
            if line.startswith("MODEL"):        # Each binding pose starts with "MODEL"
                if data == []:                  # skip 1st "MODEL"
                    continue
                df = pd.DataFrame(data,columns=pdbqtColumns)
                df[["ATOM_ID","RES_ID"]] = df[["ATOM_ID","RES_ID"]].astype(int)
                df[["X","Y","Z","OCCUPANCY","BETAFACTOR"]]=df[["X","Y","Z","OCCUPANCY","BETAFACTOR"]].astype(float)
                dockingDfList.append(df)
                data = []
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                record = [line[start:end].strip() for start, end in columsNums]
                data.append(record)
    # deal with last entry in pdbqtfile
    df = pd.DataFrame(data,columns=pdbqtColumns)
    df[["ATOM_ID","RES_ID"]] = df[["ATOM_ID","RES_ID"]].astype(int)
    df[["X","Y","Z","OCCUPANCY","BETAFACTOR"]]=df[["X","Y","Z","OCCUPANCY","BETAFACTOR"]].astype(float)
    dockingDfList.append(df)

    return dockingDfList
#########################################################################################################################
def run_vina(outDir,configFile):
    logFile = p.join(outDir,"vina_docking.log")
    with open(logFile,"a") as logFile:
        call(["vina","--config",configFile],stdout=logFile)
#########################################################################################################################
## writes a config file for a Vina docking simulation
def write_vina_config(outDir,receptorPdbqt,ligandPdbqt,boxCenter,boxSize,flexPdbqt=None,
                        exhaustiveness = 16, numModes = 10, cpus=2, energyRange = 5, seed = 42, flex=False):
    vinaConfigFile=p.join(outDir,f"vina_conf.txt")

    with open(vinaConfigFile,"w") as outFile:
        if not flex:
            outFile.write(f"receptor = {receptorPdbqt}\n")
            outFile.write(f"ligand = {ligandPdbqt}\n\n")
        else:
            outFile.write(f"receptor = {receptorPdbqt}\n")
            outFile.write(f"flex = {flexPdbqt}\n\n")
            outFile.write(f"ligand = {ligandPdbqt}\n\n")

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
            dockedPdbqt  =   p.join(outDir,f"binding_poses.pdbqt")
            outFile.write(f"out = {dockedPdbqt}\n")
        else:
            dockedPdbqt  =   p.join(outDir,f"binding_poses.pdbqt")
            outFile.write(f"out = {dockedPdbqt}\n")
        outFile.write(f"cpu = {cpus}")

        return vinaConfigFile, dockedPdbqt


#########################################################################################################################
def pocket_residues_to_alainine(protName, pdbFile, residuesToAlanine, outDir):
    residuesToAlanine = [residue for residue in residuesToAlanine]

    protDf = pdbUtils.pdb2df(pdbFile)

    pocketResidues_set = set((d['CHAIN_ID'], d['RES_ID']) for d in residuesToAlanine)

    def change_res_name(row):
        if (row['CHAIN_ID'], row['RES_ID']) in pocketResidues_set:
            return 'ALA'
        else:
            return row['RES_NAME']
    protDf['RES_NAME'] = protDf.apply(change_res_name, axis=1)

    tmpPdb = p.join(outDir, f"{protName}_tmp.pdb")
    pdbUtils.df2pdb(protDf, tmpPdb)

    firstResidue = int(protDf.iloc[0]["RES_ID"])
    protAmpal = ampal.load_pdb(tmpPdb)
    seqlength = len(protAmpal.sequences[0])
    newSequence = ''
    for i in range(firstResidue,firstResidue+seqlength):
        if i in residuesToAlanine:
            newSequence+='A'
        else:
            newSequence+=protAmpal.sequences[0][i-firstResidue]
    alaAmpal = modelling.pack_side_chains_scwrl(protAmpal,[newSequence])
    alaPdbString = alaAmpal.make_pdb(ligands=False)
    alaPdb = p.join(outDir,f"{protName}_pocketAla.pdb")
    with open(alaPdb,"w") as file:
        file.write(alaPdbString)
    return alaPdb
#########################################################################################################################
def pdb_to_pdbqt(inPdb, outDir, jobType):
    name = p.splitext(p.basename(inPdb))[0]
    outPdbqt = p.join(outDir, f"{name}.pdbqt")
    if jobType == "flex":
        obabelCommand = ["obabel", "-i","pdb", inPdb, "-o", "pdbqt", "-O", outPdbqt, "-xs"]
    elif jobType == "rigid":
        obabelCommand = ["obabel", "-i","pdb", inPdb, "-o", "pdbqt", "-O", outPdbqt, "-xr"]
    elif jobType == "ligand":
        obabelCommand = ["obabel", "-i","pdb", inPdb, "-o", "pdbqt", "-O", outPdbqt, "-xn"]
    call(obabelCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    with open(outPdbqt,"r") as f:
        fileContents = f.read()
    fileContents = fileContents.replace("NA", "N ")

    with open(outPdbqt,"w") as f:
        f.write(fileContents)


    return outPdbqt

#########################################################################################################################
def set_up_directory(protName, protDir, ligandName, ligandDir, outDir):
    runDir = p.join(outDir,f"{protName}_{ligandName}")
    os.makedirs(runDir,exist_ok=True)

    ## copy protPdb to runDir
    protPdb = p.join(protDir, f"{protName}.pdb")
    if not p.isfile(protPdb):
        print(f"Can't find pdb file:\n\t{protPdb}")
        return
    copy(protPdb,runDir)
    protPdb = p.join(runDir, f"{protName}.pdb")

    ## copy ligandPdb to runDir
    ligandPdb = p.join(ligandDir,f"{ligandName}.pdb")
    if not p.isfile(ligandPdb):
        print(f"Can't find pdb file:\n\t{ligandPdb}")
        return
    copy(ligandPdb,runDir)
    ligandPdb = p.join(runDir,f"{ligandName}.pdb")

    return protPdb, ligandPdb, runDir
#########################################################################################################################
def directed_fpocket(protName,runDir,pdbFile, targetPocketResidues):
    #print("----->\tRunning Fpocket!")
    os.chdir(runDir)
    minSphereSize = "3.0"
    maxSphereSize = "6.0"
    call(["fpocket","-f",pdbFile,"-m",minSphereSize,"-M",maxSphereSize], stdout=subprocess.PIPE)
    fpocketOutDir = p.join(runDir,f"{protName}_out","pockets")
    targetResidueCounts = {}
    ## look at all the pockets that FPocket finds, count how many target residues are in each pocket
    for file in os.listdir(fpocketOutDir):
        targetResidueCount = 0
        fileData = p.splitext(file)
        if not fileData[1] == ".pdb":
            continue
        pocketPdb = p.join(fpocketOutDir,file)
        pocketDf = pdbUtils.pdb2df(pocketPdb)
        for targetRes in targetPocketResidues:
            resData = targetRes.split(":")
            targetChain = resData[0]
            targetResId = int(resData[2])

            targetDf = pocketDf[(pocketDf["CHAIN_ID"] == targetChain) &
                                 (pocketDf["RES_ID"] == targetResId)]
            if len(targetDf) > 0:
                targetResidueCount += 1
        pocketId = fileData[0].split("_")[0]
        targetResidueCounts.update({pocketId: targetResidueCount})


    ## select the pocket with the most active site residues
    bindingPocketId = max(targetResidueCounts, key = targetResidueCounts.get)
    bindingPocketPdb = p.join(fpocketOutDir,f"{bindingPocketId}_atm.pdb")
    bindingPocketDf = pdbUtils.pdb2df(bindingPocketPdb)


    boxCenter = [bindingPocketDf["X"].mean(), bindingPocketDf["Y"].mean(),bindingPocketDf["Z"].mean()]

    pocketChains = bindingPocketDf["CHAIN_ID"].to_list()
    pocketResNames = bindingPocketDf["RES_NAME"].to_list()
    pocketResIds = bindingPocketDf["RES_ID"].to_list()

    uniqueResidues = set(zip(pocketChains, pocketResNames, pocketResIds))
    pocketResidues = [{"CHAIN_ID":chainId, "RES_NAME": resName, "RES_ID" :resId} for
                       chainId, resName, resId in uniqueResidues]
    pocketResidues = sorted(pocketResidues, key=lambda x: (x['CHAIN_ID'], x['RES_ID']))

    ## dump to yaml
    pocketResiduesYaml = p.join(runDir,f"{protName}_pocket_residues.yaml")
    with open(pocketResiduesYaml,"w") as yamlFile:
        yaml.dump(pocketResidues, yamlFile, default_flow_style=False)

    #print("----->\tFpocket Success!")

    return boxCenter, pocketResidues
#########################################################################################################################



