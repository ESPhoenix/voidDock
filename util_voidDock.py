
## import basic libraries
from itertools import chain
import os
from subprocess import call
import os.path as p
from shutil import copy, move
import pandas as pd
import json
import ampal 
import isambard.specifications as specifications 
import isambard.modelling as modelling
####################################################################################################################
def read_docking_orders(ligandOrdersCsv):
    ordersDf = pd.read_csv(ligandOrdersCsv,index_col="ID")
    ordersDf.index = ordersDf.index.astype(str)
    ordersDf = ordersDf["Ligand"]
    ordersDict = ordersDf.to_dict()

    return ordersDict
#########################################################################################################################
def get_pdb_list(protDir):
    pdbFiles =[]
    # loop through receptor PDB files
    for fileName in os.listdir(protDir):
        # Extract file name, skip if not a PDB file
        fileData = p.splitext(fileName)
        if not fileData[1] == ".pdb":
            continue
        pdbFiles.append(fileName)
    return pdbFiles
#########################################################################################################################
#########################################################################################################################
def process_vina_results(outDir,dockedPdbqt,receptorPdbqt,flex=False):
    # read output pdbqt file into a list of dataframes
    dockingDfList = read_docking_results(dockedPdbqt)
    receptorDf = pdbqt2df(receptorPdbqt)
    
    splice_docking_results(dockingDfList, receptorDf,outDir)

#########################################################################################################################
def splice_docking_results(dockingDfList, receptorDf, outDir):
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
        saveFile = p.join(finalPdbDir, f"docked_pose_{str(poseNumber)}.pdb")
        df_to_pdb(df=wholeDf,outFile=saveFile)

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
    residuesToAlanine = [int(res) for res in residuesToAlanine]

    protDf = pdb2df(pdbFile)
    firstResidue = int(protDf.iloc[0]["RES_ID"])

    protAmpal = ampal.load_pdb(pdbFile)
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
def pdb_to_pdbqt(name, pdbFile, outDir, util24Dir, mglToolsDir,jobType,flexRes=None):
    env = os.environ.copy()
    env["PYTHONPATH"] = mglToolsDir
    os.chdir(outDir)
    prepReceptorPy = p.join(util24Dir, "prepare_receptor4.py")
    prepligandPy = p.join(util24Dir,"prepare_ligand4.py")

    if jobType == "rigid":
        protPdbqt = p.join(outDir,"{}.pdbqt".format(name))
        call(["python2.7",prepReceptorPy,"-r",pdbFile,"-o",protPdbqt],env=env)
        return protPdbqt

    
    elif jobType == "ligand":
        ligandPdbqt = p.join(outDir,f"{name}.pdbqt")
        call(["python2.7",prepligandPy,"-l",pdbFile,"-o",ligandPdbqt],env=env)
        return ligandPdbqt

#########################################################################################################################
def set_up_directory(fileName,protDir,ligandDir,outDir,ordersDict):
    fileData = p.splitext(fileName)
    protPdb = p.join(protDir,fileName)
    protName = fileData[0]
    runDir = p.join(outDir,protName)
    os.makedirs(runDir,exist_ok=True)
    copy(protPdb,runDir)
    protPdb = p.join(runDir,f"{protName}.pdb")

    # identify ligand 
    ligandName = ordersDict[protName]
    ligandPdb = p.join(ligandDir,f"{ligandName}.pdb")
    copy(ligandPdb,runDir)
    ligandPdb = p.join(runDir,f"{ligandName}.pdb")

    return protName, protPdb, ligandPdb, ligandName, runDir
#########################################################################################################################
def run_fpocket(name,runDir,pdbFile):
    #print("----->\tRunning Fpocket!")
    os.chdir(runDir)
    minSphereSize = "3.0"
    maxSphereSize = "6.0"
    call(["fpocket","-f",pdbFile,"-m",minSphereSize,"-M",maxSphereSize])
    fpocketOutDir = p.join(runDir,f"{name}_out","pockets")
    ## ASSUMPTION == LARGEST POCKET IS OUR BINDING POCKET ## Not really true!
    largestPocketPdb = p.join(fpocketOutDir,"pocket1_atm.pdb")
    ## ERROR Handling
    if not p.isfile(largestPocketPdb):
        #print("--X-->\tFpocket Failed!")
        return
    largestPocketDf = pdb2df(largestPocketPdb)

    boxCenter = [largestPocketDf["X"].mean(), largestPocketDf["Y"].mean(),largestPocketDf["Z"].mean()]
    pocketResidues = largestPocketDf["RES_ID"].unique().tolist()

    #print("----->\tFpocket Success!")

    extract_pocket_info(runDir=runDir, fpocketOutDir=fpocketOutDir)    
    return boxCenter, pocketResidues

def extract_pocket_info(runDir,fpocketOutDir):
    ## LOOP THROUGH POCKET FILES ##
    pocketResidueDict = {}
    for file in os.listdir(fpocketOutDir):
        #print(file)
        ## SKIP VERT FILES ##
        if not p.splitext(file)[1] == ".pdb":
            continue
        filePath = p.join(fpocketOutDir,file)
        pocketId = file.split("_")[0]
        pocketDf = pdb2df(filePath)
        uniqueResIds = pd.unique(pocketDf["RES_ID"]).tolist()
        pocketResidueDict.update({pocketId:uniqueResIds})
    jsonDumpFile = p.join(runDir,"pocket_residues_report.json")
    with open(jsonDumpFile,"w") as jsonFile:
        json.dump(pocketResidueDict,jsonFile)

#########################################################################################################################
# read pdb files as pandas dataframes
def pdb2df(protPdb):
    columns = ['ATOM', 'ATOM_ID', 'ATOM_NAME', 'RES_NAME', 'CHAIN_ID', 'RES_ID', 'X', 'Y', 'Z', 'OCCUPANCY', 'BETAFACTOR', 'ELEMENT']
    data = []
    with open(protPdb, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_type = line[0:6].strip()
                atom_id = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21:22].strip()
                if chain_id == '':
                    chain_id = None
                res_id = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                occupancy = float(line[54:60].strip())
                temp_factor = float(line[60:66].strip())
                element = line[76:78].strip()

                data.append([atom_type, atom_id, atom_name, res_name, chain_id, res_id, x, y, z, occupancy, temp_factor, element])

    return pd.DataFrame(data, columns=columns)
##########################
# reads a pdbqt file to pandas dataframe
def pdbqt2df(pdbqtFile):
    # remove ROOT/BRANCH    
    pdbqtColumns    =   ["ATOM","ATOM_ID", "ATOM_NAME", "RES_NAME",
                    "CHAIN_ID", "RES_ID", "X", "Y", "Z", "OCCUPANCY", 
                    "BETAFACTOR","CHARGE", "ELEMENT"]
    columsNums = [(0, 6), (6, 11), (11, 17), (17, 21), (21, 22), (22, 26), 
                  (26, 38), (38, 46), (46, 54), (54, 60), (60, 70), (70, 77), (77, 79)]

    # read pdbqt file        
    data = []
    with open(pdbqtFile, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                record = [line[start:end].strip() for start, end in columsNums]
                data.append(record)
    df = pd.DataFrame(data,columns=pdbqtColumns)
    # set appropriate types for elements in dataframe
    df[["ATOM_ID","RES_ID"]] = df[["ATOM_ID","RES_ID"]].astype(int)
    df[["X","Y","Z","OCCUPANCY","BETAFACTOR"]]=df[["X","Y","Z","OCCUPANCY","BETAFACTOR"]].astype(float)
    return df
##########################
def df_to_pdb(df, outFile):
    with open(outFile,"w") as f:
        for _, row in df.iterrows():
            pdbLine = f"{row['ATOM']:<6}"
            pdbLine += f"{row['ATOM_ID']:>5}{' '*2}"
            pdbLine += f"{row['ATOM_NAME']:<4}"
            pdbLine += f"{row['RES_NAME']:<4}"
            pdbLine += f"{row['CHAIN_ID']:<1}{' '*1}"
            pdbLine += f"{row['RES_ID']:<7}"
            pdbLine += f"{row['X']:>8.3f}"
            pdbLine += f"{row['Y']:>8.3f}"
            pdbLine += f"{row['Z']:>8.3f}"
            pdbLine += f"{row['OCCUPANCY']:>6.2f}"
            pdbLine += f"{row['BETAFACTOR']:>6.2f}"
            pdbLine += "\n"
            #pdbLine += f"{row['ELEMENT']:>12}\n"
            f.write(pdbLine)
