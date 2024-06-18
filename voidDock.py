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
import config_checker
# clean code
from typing import Union
from os import PathLike


##########################################################################


def main() -> None:
    # read config file
    config: dict = config_checker.read_and_validate_config()
    outDir: Union[PathLike, str] = config["pathInfo"]["outDir"]
    ligandDir: Union[PathLike, str] = config["pathInfo"]["ligandDir"]
    dockingOrders: dict = config["dockingOrders"]

    # pre-prepare ligand pdbqt files
    gen_ligand_pdbqts(dockingOrders, ligandDir)
    # make outDir
    os.makedirs(outDir, exist_ok=True)

    cpusPerRun: int = config["cpuInfo"]["cpusPerRun"]
    parallelCpus: int = config["cpuInfo"]["totalCpuUsage"] // cpusPerRun

    if parallelCpus == 1:
        run_serial(config, dockingOrders)
    elif parallelCpus > 1:
        run_parallel(config, dockingOrders)
    collate_docked_pdbs(outDir, rmDirs = False)
##########################################################################


def run_parallel(config: dict, dockingOrders: dict) -> None:
    cpusPerRun: int = config["cpuInfo"]["cpusPerRun"]
    parallelCpus: int = config["cpuInfo"]["totalCpuUsage"] // cpusPerRun

    with mp.Pool(processes=parallelCpus) as pool:
        try:
            pool.starmap(
                docking_protocol, [
                    (config, dockingOrder) for dockingOrder in dockingOrders])
        except Exception as e:
            print(f"ERROR: {e}")
##########################################################################


def run_serial(config: dict, dockingOrders: dict) -> None:
    for dockingOrder in dockingOrders:
        docking_protocol(config, dockingOrder)


##########################################################################
def docking_protocol(config: dict, dockingOrder:dict) -> None:
    # read config
    pathInfo: Union[PathLike, str] = config["pathInfo"]
    outDir: Union[PathLike, str] = pathInfo["outDir"]
    cpusPerRun: Union[PathLike, str] = config["cpuInfo"]["cpusPerRun"]


    # set up run directory and output key variables
    protName, protPdb, ligPdbqts, runDir = set_up_directory(outDir=outDir,
                                                             pathInfo=pathInfo,
                                                               dockingOrder=dockingOrder)

    # Use fpocket to identify correct pocket, calclate box center and return
    # residues in pocket
    targetPocketResidues: list = dockingOrder["pocketResidues"]
    boxCenter, pocketResidues = directed_fpocket(protName=protName,
                                                 runDir=runDir,
                                                 pdbFile=protPdb,
                                                 targetPocketResidues=targetPocketResidues)


    if dockingOrder["mutatePocketToAla"]: 
        # Replace pocket residues with alanine
        protPdb: Union[PathLike, str] = pocket_residues_to_alainine(protName=protName,
                                            pdbFile=protPdb,
                                            residuesToAlanine=pocketResidues,
                                            dockingOrder=dockingOrder,
                                            outDir=runDir)

    flexibleDocking: bool = False
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
        rigidPdbqt: Union[PathLike, str] = pdb_to_pdbqt(inPdb=protPdb,
                                outDir=runDir,
                                jobType="rigid")

        # Write a config file for vina
        vinaConfig, dockedPdbqt = write_vina_config(outDir=runDir,
                                            receptorPdbqt=rigidPdbqt,
                                            boxCenter=boxCenter,
                                            boxSize=30,
                                            cpus=cpusPerRun,
                                            flex=False)


    # Run vina docking
    run_vina(outDir=runDir,
             configFile=vinaConfig,
             ligPdbqts=ligPdbqts)

        
    process_vina_results(outDir= runDir,
                        dockedPdbqt = dockedPdbqt,
                        receptorPdbqt = rigidPdbqt,
                        dockingOrder = dockingOrder,
                        ligandDir= config["pathInfo"]["ligandDir"]
                        )

##########################################################################
if __name__ == "__main__":
    main()
