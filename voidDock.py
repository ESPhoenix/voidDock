#########################################################################################################################
## import basic libraries
import os
import os.path as p
import sys
import argpass
import multiprocessing as mp

## inport util scripts
from util_voidDock import *
#########################################################################################################################
# get inputs
def read_inputs():
    # create an argpass parser, read config file, snip off ".py" if on the end of file
    parser = argpass.ArgumentParser()
    parser.add_argument("--config")
    args = parser.parse_args()
    configName=args.config
    configName = p.splitext(configName)[0]

    # add config to PYTHONPATH
    cwd = os.getcwd()
    configPath = p.join(cwd,configName)
    sys.path.append(configPath)
    # import config file and run input function to return variables
    try:
        config_module = __import__(configName)
        (protDir, ligandDir, outDir, mglToolsDir, util24Dir, ligandOrdersCsv) = config_module.inputs()
        return (protDir, ligandDir, outDir, mglToolsDir, util24Dir, ligandOrdersCsv)
    except ImportError:
        print(f"Error: Can't to import module '{configName}'. Make sure the input exists!")
        print("HOPE IS THE FIRST STEP ON THE ROAD TO DISAPPOINTMENT")
        exit()

#########################################################################################################################
def main():
    protDir, ligandDir, outDir, mglToolsDir, util24Dir, ligandOrdersCsv = read_inputs()

    # make outDir
    os.makedirs(outDir,exist_ok=True)    
    # read ligand orders csv file into a dictionary
    ordersDict = read_docking_orders(ligandOrdersCsv=ligandOrdersCsv)
    # get a list of receptor pdbfiles in receptor directory
    pdbFiles = get_pdb_list(protDir=protDir)
    # run docking with multiprocessing
    run_with_multiprocessing(pdbFiles,protDir, ligandDir,outDir,ordersDict,util24Dir,mglToolsDir)
    # run in serial (good for debugging. hash out otherwise!)
    #run_serial(pdbFiles,protDir, ligandDir,outDir,ordersDict,util24Dir,mglToolsDir)
#########################################################################################################################
def run_with_multiprocessing(pdbFiles,protDir, ligandDir,outDir,ordersDict,util24Dir,mglToolsDir):
    numCores = mp.cpu_count()
    with mp.Pool(processes=round(numCores / 2)) as pool:
        try:
            pool.starmap(docking_protocol, [(fileName,protDir, ligandDir,outDir,
                                             ordersDict,util24Dir,mglToolsDir) for fileName in pdbFiles])
        except Exception as e:
            print(f"ERROR: {e}")
#########################################################################################################################
def run_serial(pdbFiles,protDir, ligandDir,outDir,ordersDict,util24Dir,mglToolsDir):
    for fileName in pdbFiles:
        docking_protocol(fileName,protDir, ligandDir,outDir,ordersDict,util24Dir,mglToolsDir)


#########################################################################################################################
def docking_protocol(fileName,protDir, ligandDir,outDir,ordersDict,util24Dir,mglToolsDir):
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
                configFile = vinaConfig)
    # split docking results PDBQT file into separate PDB files
    process_vina_results(outDir = runDir,
                            receptorPdbqt = alaPdbtq,
                            dockedPdbqt = dockedPdbqt)

#########################################################################################################################
main()