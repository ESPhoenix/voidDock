##########################################################################
# import basic libraries
import os
import os.path as p
import sys
import argpass
import multiprocessing as mp
import yaml
# inport util scripts
from modules_voidDock import gnina_protocol_voidDock as gnina_protocol
from modules_voidDock import vina_protocol_voidDock as vina_protocol
from modules_voidDock import vina_prep_voidDock as vina_prep
from modules_voidDock import gnina_prep_voidDock as gnina_prep
from  modules_voidDock import config_checker as config_checker
# clean code
from typing import Union, Tuple
from os import PathLike


##########################################################################


def main() -> None:
    # read config file
    config: dict = config_checker.read_and_validate_config()
    outDir: Union[PathLike, str] = config["pathInfo"]["outDir"]
    ligandDir: Union[PathLike, str] = config["pathInfo"]["ligandDir"]
    dockingOrders: dict = config["dockingOrders"]


    # make outDir
    os.makedirs(outDir, exist_ok=True)

    cpusPerRun: int = config["cpuInfo"]["cpusPerRun"]
    parallelCpus: int = config["cpuInfo"]["totalCpuUsage"] // cpusPerRun


    ## read methodInfo and decide to use VINA or GNINA
    methodInfo: dict = config["methodInfo"]
    if methodInfo["dockingMethod"].upper() == "VINA":
        dockingProtocolFunc = vina_protocol.docking_protocol
        # pre-prepare ligand pdbqt files
        vina_prep.gen_ligand_pdbqts(dockingOrders, ligandDir)
    elif methodInfo["dockingMethod"].upper() == "GNINA":
        dockingProtocolFunc = gnina_protocol.docking_protocol
        # pre-prepare ligand sdf files
        gnina_prep.gen_ligand_sdfs(dockingOrders, ligandDir)
    
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
