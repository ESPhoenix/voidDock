##########################################################################
# import basic libraries
import os
import os.path as p
import sys
import argpass
import multiprocessing as mp
import yaml
# inport util scripts
from modules_voidDock import gnina_protocol_voidDock as gninaProtocol
from modules_voidDock import vina_protocol_voidDock as vinaProtocol
from modules_voidDock import shared_prep_voidDock as sharedPrep
from modules_voidDock import config_checker as configChecker
# clean code
from typing import Union, Tuple
from os import PathLike


##########################################################################


def main() -> None:
    # read config file
    config: dict = configChecker.read_and_validate_config()
    outDir: Union[PathLike, str] = config["pathInfo"]["outDir"]
    ligandDir: Union[PathLike, str] = config["pathInfo"]["ligandDir"]
    dockingOrders: dict = config["dockingOrders"]


    # make outDir
    os.makedirs(outDir, exist_ok=True)

    cpusPerRun: int = config["cpuInfo"]["cpusPerRun"]
    parallelCpus: int = config["cpuInfo"]["totalCpuUsage"] // cpusPerRun

    # pre-prepare ligand pdbqt files    
    sharedPrep.gen_ligand_pdbqts(dockingOrders, ligandDir)

    ## read methodInfo and decide to use VINA or GNINA
    methodInfo: dict = config["methodInfo"]
    if methodInfo["dockingMethod"].upper() == "VINA":
        dockingProtocolFunc = vinaProtocol.docking_protocol
    elif methodInfo["dockingMethod"].upper() == "GNINA":
        dockingProtocolFunc = gninaProtocol.docking_protocol
    
    ## decide to run serial or parallel
    if parallelCpus == 1:
        run_serial(config, dockingOrders, dockingProtocolFunc)
    elif parallelCpus > 1:
        run_parallel(config, dockingOrders, dockingProtocolFunc)
    # collate_docked_pdbs(outDir, rmDirs = False)
##########################################################################


def run_parallel(config: dict, dockingOrders: dict, dockingProtocolFunc: callable) -> None:
    cpusPerRun: int = config["cpuInfo"]["cpusPerRun"]
    parallelCpus: int = config["cpuInfo"]["totalCpuUsage"] // cpusPerRun

    with mp.Pool(processes=parallelCpus) as pool:
        try:
            pool.starmap(
                dockingProtocolFunc, [
                    (config, dockingOrder) for dockingOrder in dockingOrders])
        except Exception as e:
            print(f"ERROR: {e}")
##########################################################################


def run_serial(config: dict, dockingOrders: dict, dockingProtocolFunc: callable) -> None:
    for dockingOrder in dockingOrders:
        dockingProtocolFunc(config, dockingOrder)



##########################################################################
if __name__ == "__main__":
    main()
