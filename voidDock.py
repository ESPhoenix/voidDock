#########################################################################################################################
## import basic libraries
import os
import os.path as p
import sys
import argpass
import multiprocessing as mp
import yaml
## inport util scripts
from util_voidDock import *
#########################################################################################################################
# get inputs
def read_inputs():
    ## create an argpass parser, read config file, 
    parser = argpass.ArgumentParser()
    parser.add_argument("--config")
    args = parser.parse_args()

    configFile=args.config
    ## Read config.yaml into a dictionary
    with open(configFile,"r") as yamlFile:
        config = yaml.safe_load(yamlFile) 

    return config
#########################################################################################################################
def main(configFile):
    ## read config file
    config = read_inputs()
    protDir = config["dockingTargetsInfo"]["protDir"]
    ligandDir = config["dockingTargetsInfo"]["ligandDir"]
    outDir = config["dockingTargetsInfo"]["outDir"]

    cpuInfo = config["cpuInfo"]

    dockingOrders = config["dockingOrders"]
    # make outDir
    os.makedirs(outDir,exist_ok=True)    

    if cpuInfo["totalCpuUsage"] == 1:
        run_serial(protDir, ligandDir,outDir, cpuInfo,  dockingOrders)
    elif cpuInfo["totalCpuUsage"] > 1:
        run_parallel(protDir, ligandDir,outDir, cpuInfo, dockingOrders)


#########################################################################################################################
def run_parallel(protDir, ligandDir,outDir, cpuInfo,  dockingOrders):
    cpusPerRun = cpuInfo["cpusPerRun"]
    parallelCpus = cpuInfo["totalCpuUsage"] // cpusPerRun
    
    with mp.Pool(processes=parallelCpus) as pool:
        try:
            pool.starmap(docking_protocol, [(protDir, ligandDir ,outDir, cpusPerRun, dockingOrder)
                                             for dockingOrder in dockingOrders])
        except Exception as e:
            print(f"ERROR: {e}")
#########################################################################################################################
def run_serial(protDir, ligandDir,outDir, cpuInfo,  dockingOrders):
    cpusPerRun = cpuInfo["cpusPerRun"]
    for dockingOrder in dockingOrders:
        docking_protocol(protDir, ligandDir ,outDir, cpusPerRun, dockingOrder)


#########################################################################################################################
def docking_protocol(protDir, ligandDir ,outDir, cpusPerRun, dockingOrder):
    ## read order for this docking run from dockingOrder
    protName = dockingOrder["protein"]
    ligandName = dockingOrder["ligand"]


    # set up run directory and output key variables
    protPdb, ligandPdb, runDir = set_up_directory(protDir = protDir, protName=protName,
                                                    ligandDir = ligandDir, ligandName = ligandName,
                                                    outDir=outDir)
    # Use fpocket to identify correct pocket, calclate box center and return residues in pocket
    targetPocketResidues = dockingOrder["pocketResidues"]
    boxCenter, pocketResidues       =   directed_fpocket(protName=protName,
                                                        runDir=runDir,
                                                        pdbFile=protPdb,
                                                        targetPocketResidues = targetPocketResidues)
    # Replace pocket residues with alanine
    alaPdb                          =    pocket_residues_to_alainine(protName=protName,
                                                                        pdbFile=protPdb,
                                                                        residuesToAlanine = pocketResidues,
                                                                        outDir=runDir)
    # Convert alanine PDB to PDBQT
    alaPdbtq                        =   pdb_to_pdbqt(inPdb=alaPdb,
                                                        outDir = runDir,
                                                        jobType = "rigid")
    # Convert ligand PDB to PDBQT
    ligandPdbqt                     =   pdb_to_pdbqt(inPdb=ligandPdb,
                                                        outDir = runDir,
                                                        jobType = "ligand")
    # Write a config file for vina
    vinaConfig, dockedPdbqt         =   write_vina_config(outDir = runDir,
                                                            receptorPdbqt = alaPdbtq,
                                                            ligandPdbqt = ligandPdbqt,
                                                            boxCenter = boxCenter,
                                                            boxSize = 30,
                                                            cpus=cpusPerRun)
    # Run vina docking
    run_vina(outDir = runDir,
                configFile = vinaConfig)
    # split docking results PDBQT file into separate PDB files
    process_vina_results(protName= protName,
                         ligandName= ligandName,
                        outDir = runDir,
                        receptorPdbqt = alaPdbtq,
                        dockedPdbqt = dockedPdbqt)

#########################################################################################################################
if __name__ == "__main__":
    main(configFile=None)
