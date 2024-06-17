from typing import Tuple, Union, 
import argpass
import yaml
import os
from os import PathLike
from os import path as p


def check_dict_for_key(info: dict, key: any) -> Tuple[bool, any]:
    if key in info:
        if info[key] == False:
            return True, False
        else:
            return True, info[key]
    return False, None

def validate_path(argName, argPath):
    if  not isinstance(argPath, Union[os.PathLike, str]) :
        raise TypeError(f"The config argument {argName} = {argPath} is not a PathLike.")
    # Check if the path exists
    if not p.exists(argPath):
        raise FileNotFoundError(f"The config argument {argName} = {argPath} does not exist.")


def check_pathInfo(config: dict) -> None:
    isPathInfo, pathInfo = check_dict_for_key(config, "pathInfo")
    if not isPathInfo:
        print("No pathInfo found in config file")
        exit(1)
    for argName in ["protDir", "ligandDir", "outDir"]:
        isArg, argValue = check_dict_for_key(pathInfo, argName)
        if not isArg:
            print(f"Argument {argName} not found in pathInfo")
        validate_path(argName, argValue)


def check_cpuInfo(config: dict) -> None:
    isCpuInfo, cpuInfo = check_dict_for_key(config, "cpuInfo")
    if not isCpuInfo:
        print("No cpuInfo found in config file")
        exit(1)
    for argName in ["totalCpuUsage", "cpusPerRun"]:
        isArg, argValue = check_dict_for_key(cpuInfo, argName)
        if not isArg:
            print(f"Argument {argName} not found in cpuInfo")
            exit(1)
        if not isinstance(argValue, int):
            raise TypeError(f"The config argument {argName} = {argValue} is not a an int type.")

def check_dockingOrders(config: dict) -> None:
    ## we know pathInfo is all good, so we can use it here
    pathInfo = config["pathInfo"]
    protDir = pathInfo["protDir"]
    ligandDir = pathInfo["ligandDir"]


    isDockingOrders, dockingOrders = check_dict_for_key(config, "dockingOrders")
    if not isDockingOrders:
        print("No dockingOrders found in config file")
        exit(1)
    if len(dockingOrders) == 0:
        print("No entries in dockingOrders argument in config")
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