##########################################################################
# import basic libraries
import os
import os.path as p
import sys
import argpass
import multiprocessing as mp
import yaml
# inport util scripts
from util_voidDock import *
##########################################################################
# get inputs


def read_inputs():
    # create an argpass parser, read config file,
    parser = argpass.ArgumentParser()
    parser.add_argument("--config")
    args = parser.parse_args()

    configFile = args.config
    # Read config.yaml into a dictionary
    with open(configFile, "r") as yamlFile:
        config = yaml.safe_load(yamlFile)

    return config
##########################################################################


def main(configFile):
    # read config file
    config = read_inputs()
    outDir = config["pathInfo"]["outDir"]
    ligandDir = config["pathInfo"]["ligandDir"]
    dockingOrders = config["dockingOrders"]

    # pre-prepare ligand pdbqt files
    gen_ligand_pdbqts(dockingOrders, ligandDir)
    # make outDir
    os.makedirs(outDir, exist_ok=True)

    cpusPerRun = config["cpuInfo"]["cpusPerRun"]
    parallelCpus = config["cpuInfo"]["totalCpuUsage"] // cpusPerRun

    if parallelCpus == 1:
        run_serial(config, dockingOrders)
    elif parallelCpus > 1:
        run_parallel(config, dockingOrders)
    # collate_docked_pdbs(outDir) ##TODO UNHASH!
##########################################################################


def run_parallel(config, dockingOrders):
    cpusPerRun = config["cpuInfo"]["cpusPerRun"]
    parallelCpus = config["cpuInfo"]["totalCpuUsage"] // cpusPerRun

    with mp.Pool(processes=parallelCpus) as pool:
        try:
            pool.starmap(
                docking_protocol, [
                    (config, dockingOrder) for dockingOrder in dockingOrders])
        except Exception as e:
            print(f"ERROR: {e}")
##########################################################################


def run_serial(config, dockingOrders):
    for dockingOrder in dockingOrders:
        docking_protocol(config, dockingOrder)


##########################################################################
def docking_protocol(config, dockingOrder):
    # read config
    pathInfo = config["pathInfo"]
    outDir = pathInfo["outDir"]
    cpusPerRun = config["cpuInfo"]["cpusPerRun"]

    # set up run directory and output key variables
    protName, protPdb, ligPdbqts, runDir = set_up_directory(
        outDir=outDir, pathInfo=pathInfo, dockingOrder=dockingOrder)

    # Use fpocket to identify correct pocket, calclate box center and return
    # residues in pocket
    targetPocketResidues = dockingOrder["pocketResidues"]
    boxCenter, pocketResidues = directed_fpocket(protName=protName,
                                                 runDir=runDir,
                                                 pdbFile=protPdb,
                                                 targetPocketResidues=targetPocketResidues)


    if dockingOrder["mutatePocketToAla"]: 
        # Replace pocket residues with alanine
        protPdb = pocket_residues_to_alainine(protName=protName,
                                            pdbFile=protPdb,
                                            residuesToAlanine=pocketResidues,
                                            dockingOrder=dockingOrder,
                                            outDir=runDir)

    flexibleDocking = False
    if "flexibleResidues" in dockingOrder:
        if len(dockingOrder["flexibleResidues"]) > 0:
            flexibleDocking  = True
            rigidPdbqt, flexPdbqt = gen_flex_pdbqts(protPdb = protPdb,
                            flexibleResidues = dockingOrder["flexibleResidues"],
                            outDir = runDir)

            # Write a config file for vina
            vinaConfig, dockedPdbqt = write_vina_config(outDir=runDir,
                                                receptorPdbqt=rigidPdbqt,
                                                flexPdbqt = flexPdbqt,
                                                boxCenter=boxCenter,
                                                boxSize=30,
                                                cpus=cpusPerRun,
                                                flex=True)
                                    
    if not flexibleDocking:
        # Convert alanine PDB to PDBQT
        protPdbqt = pdb_to_pdbqt(inPdb=protPdb,
                                outDir=runDir,
                                jobType="rigid")

        # Write a config file for vina
        vinaConfig, dockedPdbqt = write_vina_config(outDir=runDir,
                                            receptorPdbqt=protPdbqt,
                                            boxCenter=boxCenter,
                                            boxSize=30,
                                            cpus=cpusPerRun,
                                            flex=False)


    # Run vina docking
    run_vina(outDir=runDir,
             configFile=vinaConfig,
             ligPdbqts=ligPdbqts)

    # split docking results PDBQT file into separate PDB files
    process_vina_results(dockingOrder=dockingOrder,
                         outDir=runDir,
                         receptorPdb=protPdb,
                         dockedPdbqt=dockedPdbqt)


##########################################################################
if __name__ == "__main__":
    main(configFile=None)
