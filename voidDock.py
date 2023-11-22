#       GOALS
#
#       Loop over all PDB files in protDir
#           --> find binding pocket with fpocket                     [x]
#           ---> convert ligand and          [ ]
#           ----> create binding poses with CAVER DOCK              [ ]
#
#
#
#
#########################################################################################################################
## import basic libraries
from itertools import chain
import os
from subprocess import call
import os.path as p
from shutil import copy, move
import pandas as pd
import sys
import ampal 
import isambard.specifications as specifications 
import isambard.modelling as modelling

import multiprocessing as mp
#########################################################################################################################
def inputs():
    protDir = "/home/eugene/DOCK_test/input_pdbs"
    ligandDir = "/home/eugene/DOCK_test/ligands"
    outDir = "/home/eugene/DOCK_test/outputs"
    mglToolsDir = "/home/eugene/bin/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs"
    util24Dir = "/home/eugene/bin/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24"
    vinaExe = "/home/eugene/bin/vina/vina_1.2.5_linux_x86_64"
    ligandOrdersCsv = "/home/eugene/DOCK_test/ligand_orders.csv"

    return protDir, ligandDir, outDir, mglToolsDir, util24Dir, vinaExe, ligandOrdersCsv
#########################################################################################################################
#########################################################################################################################
def main():
    protDir, ligandDir, outDir, mglToolsDir, util24Dir, vinaExe, ligandOrdersCsv = inputs()

    # make outDir
    ifnotmkdir(outDir)    
    # read ligand orders csv file into a dictionary
    ordersDf = pd.read_csv(ligandOrdersCsv,index_col="ID")
    ordersDf = ordersDf["Cofactor"]
    ordersDict = ordersDf.to_dict()

    pdbFiles =[]
    # loop through receptor PDB files
    for fileName in os.listdir(protDir):
        # Extract file name, skip if not a PDB file
        fileData = p.splitext(fileName)
        if not fileData[1] == ".pdb":
            continue
        pdbFiles.append(fileName)

    num_cores = mp.cpu_count()
    pool = mp.Pool(processes=round(num_cores / 2))

    for fileName in pdbFiles:
        pool.apply_async(docking_protocol, (fileName,protDir, ligandDir,outDir,ordersDict,util24Dir,mglToolsDir,vinaExe))

    # Close the pool to release resources
    pool.close()
    pool.join()

#########################################################################################################################
#########################################################################################################################

def docking_protocol(fileName,protDir, ligandDir,outDir,ordersDict,util24Dir,mglToolsDir,vinaExe):
    # set up run directory and output key variables
    protName, protPdb, ligandPdb, ligandName, runDir = set_up_directory(fileName=fileName,
                                                                            protDir=protDir,
                                                                            ligandDir = ligandDir,
                                                                            outDir=outDir,
                                                                            ordersDict=ordersDict)  
    # Use fpocket to identify largest pocket, return center of pocket as [X,Y,Z] coords and p=pocket residues
    boxCenter, pocketResidues       =   run_fpocket(name=protName,
                                                        runDir=runDir,
                                                        pdbFile=protPdb)
    # Replace pocket residues with alanine
    alaPdb                          =    pocket_residues_to_alainine(protName=protName,
                                                                        pdbFile=protPdb,
                                                                        residuesToAlanine = pocketResidues,
                                                                        outDir=runDir)
    # Convert alanine PDB to PDBQT
    alaPdbtq                        =   pdb_to_pdbqt(name = protName,
                                                        pdbFile=alaPdb,
                                                        outDir = runDir,
                                                        util24Dir = util24Dir,
                                                        mglToolsDir  = mglToolsDir,
                                                        jobType = "rigid")
    # Convert ligand PDB to PDBQT
    ligandPdbqt                     =   pdb_to_pdbqt(name = ligandName,
                                                        pdbFile=ligandPdb,
                                                        outDir = runDir,
                                                        util24Dir = util24Dir,
                                                        mglToolsDir  = mglToolsDir,
                                                        jobType = "ligand")
    # Write a config file for vina
    vinaConfig, dockedPdbqt         =   write_vina_config(outDir = runDir,
                                                            receptorPdbqt = alaPdbtq,
                                                            ligandPdbqt = ligandPdbqt,
                                                            boxCenter = boxCenter,
                                                            boxSize = 30)
    # Run vina docking
    run_vina(outDir = runDir,
                configFile = vinaConfig,
                vinaExe = vinaExe)
    # split docking results PDBQT file into separate PDB files
    process_vina_results(outDir = runDir,
                            receptorPdbqt = alaPdbtq,
                            dockedPdbqt = dockedPdbqt)


#########################################################################################################################
#########################################################################################################################
def process_vina_results(outDir,dockedPdbqt,receptorPdbqt):
    # remove ROOT/BRANCH
    pdbqtColumns    =   ["ATOM","ATOM_ID", "ATOM_NAME", "RES_NAME",
                    "CHAIN_ID", "RES_ID", "X", "Y", "Z", "OCCUPANCY", 
                    "BETAFACTOR","CHARGE", "ELEMENT"]
    columsNums = [(0, 6), (6, 11), (11, 17), (17, 21), (21, 22), (22, 26), 
                  (26, 38), (38, 46), (46, 54), (54, 60), (60, 70), (70, 77), (77, 79)]
    # read pdbqt file into multiple dataframes
    dockingDfList =[]
    data = []

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



    receptorDf = pdbqt2df(receptorPdbqt)
    finalPdbDir = p.join(outDir,"final_docked_pdbs")
    ifnotmkdir(finalPdbDir)
    ligandResidueId = str(int(receptorDf.iloc[-1]["RES_ID"])+1)
    lastChainIdInProt = receptorDf.iloc[-1]["CHAIN_ID"]
    ligandChainId = chr((ord(lastChainIdInProt) - ord('A') + 1) % 26 + ord('A'))

    for poseNumber, dockedDf in zip(range(1,len(dockingDfList)+1),dockingDfList):
        dockedDf["RES_ID"] = ligandResidueId
        dockedDf["ATOM"] = "ATOM"
        dockedDf["CHAIN_ID"] = ligandChainId
        newDf = pd.concat([receptorDf,dockedDf],axis=0)
        newDf["ATOM_ID"] = range(1,len(newDf)+1)

        writePdbFile = p.join(finalPdbDir, f"docked_pose_{str(poseNumber)}.pdb")
        df_to_pdb(df = newDf,outFile=writePdbFile)

#########################################################################################################################
def run_vina(outDir,configFile,vinaExe):
    logFile = p.join(outDir,"vina_docking.log")
    with open(logFile,"a") as logFile:
        call([vinaExe,"--config",configFile],stdout=logFile)
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
def select_flexible_residues(protName,protPdb,flexResList):
    # load protein PDB as a DF
    protDf = pdb2df(protPdb)
    # reduce DF down to unique residues in flexResList (created from binding pocket)
    flexDf = protDf[protDf["RES_ID"].isin(flexResList)]
    flexDf.drop_duplicates(subset=["RES_ID"], inplace=True)
    # remove static residues from flexible residues (no dihedrals to rotate anyway!)
    staticResidues = ["ALA","PRO","GLY"]
    flexDf = flexDf[~flexDf["RES_NAME"].isin(staticResidues)]
    # get a list of unique CHAIN IDs
    uniqueChains = flexDf.drop_duplicates(subset=["CHAIN_ID"])["CHAIN_ID"].tolist()
    # generate MGLTools compatable flexible residues
    flexResidues=[]
    for chain in uniqueChains:
        chainBlock=[]
        for _, row in flexDf.iterrows():
            if row["CHAIN_ID"]==chain:
                resName = row["RES_NAME"]
                resSeq = row["RES_ID"]
                selection=f"{resName}{resSeq}"
                chainBlock.append(selection)
        chainBlock="_".join(chainBlock)
        chainBlock=f"{protName}:{chain}:{chainBlock}"
        flexResidues.append(chainBlock)
    
    flexResidues=",".join(flexResidues)
    return flexResidues            

#########################################################################################################################

def pdb_to_pdbqt(name, pdbFile, outDir, util24Dir, mglToolsDir,jobType,flexRes=None):
    env = os.environ.copy()
    env["PYTHONPATH"] = mglToolsDir
    os.chdir(outDir)
    prepReceptorPy = p.join(util24Dir, "prepare_receptor4.py")
    prepligandPy = p.join(util24Dir,"prepare_ligand4.py")
    prepFlexReceptorPy = p.join(util24Dir,"prepare_flexreceptor4.py")

    if jobType == "rigid":
        protPdbqt = p.join(outDir,"{}.pdbqt".format(name))
        call(["python2",prepReceptorPy,"-r",pdbFile,"-o",protPdbqt],env=env)
        return protPdbqt

    elif jobType == "flex":
        if flexRes == None:
            print(f"--X-->\tNo Flexible residues supplied for {name}")
            return
        rigidPdbqt = p.join(outDir,f"{name}_rigid.pdbqt")
        flexPdbqt = p.join(outDir,f"{name}_flex.pdbqt")
        call(["python2",prepReceptorPy,"-r",pdbFile,"-s",flexRes,"-g",rigidPdbqt,"-x",flexPdbqt],env=env)
        return rigidPdbqt, flexPdbqt
    
    elif jobType == "ligand":
        ligandPdbqt = p.join(outDir,f"{name}.pdbqt")
        call(["python2",prepligandPy,"-l",pdbFile,"-o",ligandPdbqt],env=env)
        return ligandPdbqt

#########################################################################################################################
def set_up_directory(fileName,protDir,ligandDir,outDir,ordersDict):
    fileData = p.splitext(fileName)
    protPdb = p.join(protDir,fileName)
    protName = fileData[0]
    runDir = ifnotmkdir(p.join(outDir,protName))
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
    print("----->\tRunning Fpocket!")
    os.chdir(runDir)
    minSphereSize = "3.0"
    maxSphereSize = "6.0"
    call(["fpocket","-f",pdbFile,"-m",minSphereSize,"-M",maxSphereSize])
    fpocketOutDir = p.join(runDir,f"{name}_out","pockets")
    ## ASSUMPTION == LARGEST POCKET IS OUR BINDING POCKET ## Not really true!
    largestPocketPdb = p.join(fpocketOutDir,"pocket1_atm.pdb")
    ## ERROR Handling
    if not p.isfile(largestPocketPdb):
        print("--X-->\tFpocket Failed!")
        return
    largestPocketDf = pdb2df(largestPocketPdb)

    boxCenter = [largestPocketDf["X"].mean(), largestPocketDf["Y"].mean(),largestPocketDf["Z"].mean()]
    pocketResidues = largestPocketDf["RES_ID"].unique().tolist()

    print("----->\tFpocket Success!")
    return boxCenter, pocketResidues




#########################################################################################################################
#########################################################################################################################
# basic dir creation
def ifnotmkdir(dir):
    if not p.isdir(dir):
        os.mkdir(dir)
    return dir
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

#########################################################################################################################
main()