from typing import Tuple, Union
import argpass
import yaml
import os
from os import PathLike
from os import path as p
import multiprocessing as mp

#########################################################################
def read_and_validate_config() -> dict:
    """
    Main script for config checker:
    1. Accepts --config flag arg with argpass, reads yaml file into a dict
    2. Checks paths in pathInfo to see if they are real
    3. Checks inputs of cpuInfo to see if they are correct
    4. Checks each dockingOrder in dockingOrders to see if they are correct

    Returns:
    - config (dict)
    """

    print("checking config file...")
    config: dict = read_input_yaml()
    check_pathInfo(config)
    check_cpuInfo(config)
    check_dockingOrders(config)
    print("... config file is correct")
    return config

#########################################################################
def read_input_yaml() -> dict:
    """
    Reads a YAML file using the "--config" flag with argpass
    Reads YAML file into a dict

    Returns:
    - config (dict)
    """
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
        raise FileNotFoundError(f"config file {configFile} not found")
    except yaml.YAMLError as exc:
        raise yaml.YAMLError("Error parsing YAML file:", exc)



#########################################################################
def check_pathInfo(config: dict) -> None:
    """
    Checks for pathInfo entry in config
    Checks paths in pathInfo to see if they are real
    Don't check outDir, this will be made automatically
    """
    ## check if pathInfo in config
    pathInfo, = check_info_for_args(config, "config", ["pathInfo"], optional=False)
    # check for required args in pathInfo
    protDir, ligandDir, outDir = check_info_for_args(pathInfo, "pathInfo", ["protDir", "ligandDir", "outDir"], optional=False)
    ## make sure paths exist
    for argValue, argName in zip([protDir, ligandDir], ["protDir", "ligandDir"]):
        validate_path(argName, argValue)
#########################################################################

def check_cpuInfo(config: dict) -> None:
    """
    Checks for cpuInfo in config
    Checks that entries in cpuInfo are int types
    Makes sure your computer has enough CPUs for the simulation
    """

    ## check for cpuInfo in config
    cpuInfo, = check_info_for_args(config, "config", ["cpuInfo"], optional= False)
    ## check for required args in cpuInfo
    totalCpuUseage, cpusPerRun = check_info_for_args(cpuInfo, "cpuInfo", ["totalCpuUsage", "cpusPerRun"], optional= False)
    ## check that cpuInfo arguments are int values and are positive
    for argValue, argName in zip([totalCpuUseage, cpusPerRun], ["totalCpuUseage", "cpusPerRun"]):
        if not isinstance(argValue, int):
            raise TypeError(f"The config argument {argName} = {argValue} is not a an int type.")
        if argValue < 1:
            raise ValueError(f"The config argument {argName} = {argValue} must be a int greater than 1")
    ## check that your computer has enough CPUs
    if totalCpuUseage > mp.cpu_count():
        raise ValueError("totalCpuUseage argument exceeds your computers number of cores")



#########################################################################
def check_dockingOrders(config: dict) -> None:
    """
    Checks for dockingOrders in config, makes sure it's a list and not empty
    For each entry in dockingOrders:
    - checks for required args
    - ensures that input files exist in the expected directories
    - checks for optional args
    - ensures that residue dicts in optional args are formatted correctly

    NOTE that this checks basic formatting - if you have specified the
    wrong residue with the correct format,  this will not ne caught.
    The docking simulation prep will fail at some point later on!
    """

    ## we know pathInfo is all good, so we can use it here
    pathInfo = config["pathInfo"]
    protDir = pathInfo["protDir"]
    ligandDir = pathInfo["ligandDir"]
    ## check for dockingOrders in config
    dockingOrders, = check_info_for_args(config, "config", ["dockingOrders"], optional= False)
    ## ensure that dockingOrders is a non-zero-length list
    if not isinstance(dockingOrders, list):
        raise TypeError("dockingOrders must be a list containing dicts")
    if len(dockingOrders) == 0:
        raise ValueError("dockingOrders must have at least one entry")
    ## look through each entry in dockingOrders
    for dockingOrder in dockingOrders:
        ## check for required args in dockingOrder
        protName, ligandNames =   check_info_for_args(dockingOrder, "dockingOrder", ["protein", "ligands"], optional= False)
        ## check to see if protein and ligand files exist in expected directories
        protPdb: Union[PathLike, str] = p.join(protDir,f"{protName}.pdb")
        if not p.isfile(protPdb):
            raise FileNotFoundError(f"Protein PDB {protPdb} does not exist")
        ligandNames: list = dockingOrder["ligands"]
        for ligandName in ligandNames:
            ligandPdb: Union[PathLike, str] = p.join(ligandDir, f"{ligandName}.pdb")
            if not p.isfile(ligandPdb):
                raise FileNotFoundError(f"Ligand PDB {ligandPdb} does not exist")
        ## check for optional args in dockingOrder
        pocketResidues, keepResidues, flexibleResidues = check_info_for_args(dockingOrder, "dockingOrder", ["pocketResidues", "keepResidues", "flexibleResidues"], optional=True)
        ## ensure that these residue-specifying dicts are formatted correctly
        for residues, dictName in zip([pocketResidues, keepResidues, flexibleResidues],["pocketResidues", "keepResidues", "flexibleResidues"]):
            if not residues: continue
            validate_residue_dict(residues, dictName)
#########################################################################
def check_info_for_args(info: dict, infoName: str,  argNames: list, optional: bool) -> list:
    """
    Simple check to see if list of keys is in a dict
    If optional is set to "False", raises a KeyError
    If all is as expected returns a list of values to be unpacked
    """
    ## init empty list to append to 
    unpackedDicts: list = []
    for  argName in argNames:
        isArg, argValue = check_dict_for_key(info, argName)
        ## if a required arg is not in the dict, raise a KeyError
        if not isArg and not optional:
            raise KeyError(f"Argument {argName} not found in {infoName}")
        unpackedDicts.append(argValue)
    return unpackedDicts
#########################################################################
def check_dict_for_key(info: dict, key: any) -> Tuple[bool, any]:
    """
    Checks to see if a key is in a dict
    If it is, return "True" and the associated value 
    """
    if key in info:
        if info[key] == False:
            return True, False
        else:
            return True, info[key]
    return False, False
#########################################################################
def validate_path(argName: str, argPath: Union[PathLike, str]) -> None:
    """
    Check to see if a path variable is indeed the correct type
    Check to see if the path exists
    """
    if  not isinstance(argPath, (os.PathLike, str)) :
        raise TypeError(f"The config argument {argName} = {argPath} is not a PathLike.")
    # Check if the path exists
    if not p.exists(argPath):
        raise FileNotFoundError(f"The config argument {argName} = {argPath} does not exist.")
#########################################################################
def validate_residue_dict(info: dict, infoName: str) -> None:
    """
    Ensures that the residue-specifying dicts are correctly formatted, eg:
    {"CHAIN_ID": "B", "RES_NAME: "GLY", "RES_ID": 123}
    """
    aminoAcidThreeLetterNames = init_amino_acid_list()
    for residue in info:
        chainId, resName, resId = check_info_for_args(residue, "residue", ["CHAIN_ID", "RES_NAME", "RES_ID"], optional=False)
        if not isinstance(chainId, str):
            raise TypeError(f"CHAIN_ID in {infoName} must be a string")
        if not isinstance(resName, str):
            raise TypeError(f"RES_NAME in {infoName} must be a string")
        if not resName in aminoAcidThreeLetterNames:
            raise ValueError(f"RES_NAME in {infoName} must be a canonical amino acid three-letter, ALL-CAPS code")
        if not isinstance(resId, int):
            raise TypeError(f"RES_ID in {infoName} must be a int")
        if resId < 1:
            raise ValueError(f"RES_ID in {infoName} must be a positive int")
#########################################################################
def init_amino_acid_list() -> list:
    """
    Creates list of the three-letter codes for the 20 amino acids
    TODO: add method for dealing with non-naturals
    """
    aminoAcidThreeLetterNames: list = [
        "ALA",  # Alanine
        "ARG",  # Arginine
        "ASN",  # Asparagine
        "ASP",  # Aspartic acid
        "CYS",  # Cysteine
        "GLN",  # Glutamine
        "GLU",  # Glutamic acid
        "GLY",  # Glycine
        "HIS",  # Histidine
        "ILE",  # Isoleucine
        "LEU",  # Leucine
        "LYS",  # Lysine
        "MET",  # Methionine
        "PHE",  # Phenylalanine
        "PRO",  # Proline
        "SER",  # Serine
        "THR",  # Threonine
        "TRP",  # Tryptophan
        "TYR",  # Tyrosine
        "VAL"   # Valine
    ]
    return aminoAcidThreeLetterNames
