from typing import Tuple, Union
import argpass
import yaml
import os
from os import PathLike
from os import path as p
import multiprocessing as mp

#########################################################################
def read_and_validate_config() -> dict:

    print("checking config file...")
    config: dict = read_input_yaml()
    check_pathInfo(config)
    check_cpuInfo(config)
    check_dockingOrders(config)
    print("... config file is correct")
    return config

#########################################################################
def read_input_yaml() -> dict:
    # create an argpass parser, read config file,
    parser = argpass.ArgumentParser()
    parser.add_argument("--config")
    args = parser.parse_args()

    configFile: Union[PathLike, str] = args.config
    # Read config.yaml into a dictionary
    try:
        with open(configFile, "r") as yamlFile:
            config: dict = yaml.safe_load(yamlFile)
            return config
    except FileNotFoundError:
        print("Configuration file not found.")
        exit(1)
    except yaml.YAMLError as exc:
        print("Error parsing YAML file:", exc)
        exit(1)


#########################################################################
def check_pathInfo(config: dict) -> None:
    pathInfo, = check_info_for_manditory_args(config, "config", ["pathInfo"])
    protDir, ligandDir, outDir = check_info_for_manditory_args(pathInfo, "pathInfo", ["protDir", "ligandDir", "outDir"])
    for argValue, argName in zip([protDir, ligandDir, outDir], ["protDir", "ligandDir", "outDir"]):
        validate_path(argName, argValue)
#########################################################################

def check_cpuInfo(config: dict) -> None:
    cpuInfo, = check_info_for_manditory_args(config,
                                             "config",
                                                 ["cpuInfo"])

    totalCpuUseage, cpusPerRun = check_info_for_manditory_args(cpuInfo,
                                                                "cpuInfo",
                                                                  ["totalCpuUsage", "cpusPerRun"])
    ## check that cpuInfo arguments are int values and are positive
    for argValue, argName in zip([totalCpuUseage, cpusPerRun], ["totalCpuUseage", "cpusPerRun"]):
        if not isinstance(argValue, int):
            raise TypeError(f"The config argument {argName} = {argValue} is not a an int type.")
        if argValue < 1:
            raise ValueError(f"The config argument {argName} = {argValue} must be a int greater than 1")
    if totalCpuUseage > mp.cpu_count():
        raise ValueError("totalCpuUseage argument exceeds your computers number of cores")



#########################################################################
def check_dockingOrders(config: dict) -> None:
    ## we know pathInfo is all good, so we can use it here
    pathInfo = config["pathInfo"]
    protDir = pathInfo["protDir"]
    ligandDir = pathInfo["ligandDir"]
    dockingOrders, = check_info_for_manditory_args(config, "config", ["dockingOrders"])

    if len(dockingOrders) == 0:
        print("No entries in dockingOrders argument in config")
        exit(1)
    for dockingOrder in dockingOrders:
        ## mandatory args
        for  argName in ["protein", "ligands", "pocketResidues"]:
            isArg, argValue = check_dict_for_key(dockingOrder, argName)
            if not isArg:
                print(f"Argument {argName} not found in cpuInfo")
                exit(1)
        ## check to see if protein and ligand files have pdb files where you would expect them
        protName = dockingOrder["protein"] 
        protPdb = p.join(protDir,f"{protName}.pdb")
        if not p.isfile(protPdb):
            raise FileNotFoundError(f"Protein PDB {protPdb} does not exist")
        ligandNames = dockingOrder["ligands"]
        for ligandName in ligandNames:
            ligandPdb = p.join(ligandDir, f"{ligandName}.pdb")
            if not p.isfile(ligandPdb):
                raise FileNotFoundError(f"Ligand PDB {ligandPdb} does not exist")

#########################################################################
def check_info_for_manditory_args(info, infoName,  argNames):
    unpackedDicts = []
    for  argName in argNames:
        isArg, argValue = check_dict_for_key(info, argName)
        if not isArg:
            raise KeyError(f"Argument {argName} not found in {infoName}")
        unpackedDicts.append(argValue)
    return unpackedDicts
#########################################################################
def check_dict_for_key(info: dict, key: any) -> Tuple[bool, any]:
    if key in info:
        if info[key] == False:
            return True, False
        else:
            return True, info[key]
    return False, None
#########################################################################
def validate_path(argName, argPath):
    if  not isinstance(argPath, (os.PathLike, str)) :
        raise TypeError(f"The config argument {argName} = {argPath} is not a PathLike.")
    # Check if the path exists
    if not p.exists(argPath):
        raise FileNotFoundError(f"The config argument {argName} = {argPath} does not exist.")
#########################################################################
def validate_residue_dicts(info):
    for residue in info:
        chainId, resName, resId = check_info_for_manditory_args(residue, "residue", ["CHAIN_ID", "RES_NAME", "RES_ID"])
#########################################################################