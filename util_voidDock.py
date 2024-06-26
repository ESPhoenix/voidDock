
# import basic libraries
import os
import subprocess
from subprocess import call, run
import os.path as p
from shutil import copy, move, rmtree
import pandas as pd
import yaml
## WellsWood Lab libraries
import ampal
import isambard.modelling as modelling
from pdbUtils import pdbUtils
## clean code
from typing import Tuple, Union, Literal
from os import PathLike
import numpy as np
##########################################################################

def gen_ligand_pdbqts(dockingOrders: dict, ligandDir: Union[PathLike, str]) -> None:
    '''
    Before running docking, generate pdbqt files for all ligands
    This saves us from generating a ligand pdbqt per docking run
    '''
    ## loop through ligand pdb files in the docking orders
    ## append these files to a list
    ## get unique entries
    allLigands: list = []
    for dockingOrder in dockingOrders:
        ligands: list = dockingOrder["ligands"]
        for ligand in ligands:
            allLigands.append(ligand)
    allLigands: list = list(set(allLigands))
    ## look in ligandsDir for these ligands
    ## call pdb_to_pdbqt to convert to pdbqt file
    for ligand in allLigands:
        ligPdb: Union[PathLike, str] = p.join(ligandDir, f"{ligand}.pdb")
        if not p.isfile(ligPdb):
            print(f"{ligPdb} not found, skipping...")
            continue
        pdb_to_pdbqt(ligPdb, ligandDir, jobType="ligand")
##########################################################################


def collate_docked_pdbs(outDir: Union[PathLike, str], rmDirs: bool = True) -> None:
    '''
    After a docking simulation has been run, get:
     - output pdb files
     - pocket_residues.yaml (useful for further pipelines)
     Move these files to a new dir
     Remove the rest of the outputs to save space
    '''
    ## loop through run directory and copy required files to final_docked_pdbs
    for dir in os.listdir(outDir):
        runDir: Union[PathLike, str]= p.join(outDir, dir)
        finalDockedPdbDir: Union[PathLike, str]= p.join(runDir, "final_docked_pdbs")
        if not p.isdir(finalDockedPdbDir):
            continue
        for file in os.listdir(finalDockedPdbDir):
            if not p.splitext(file)[1] == ".pdb":
                continue
            pdbFile: Union[PathLike, str]= p.join(finalDockedPdbDir, file)
            pdbDest: Union[PathLike, str]= p.join(outDir, file)
            copy(pdbFile, pdbDest)
        for file in os.listdir(runDir):
            if file.endswith("_pocket_residues.yaml"):
                pocketResiduesYaml: Union[PathLike, str]= p.join(runDir, file)
                pocketResiduesDest: Union[PathLike, str]= p.join(outDir, file)
                copy(pocketResiduesYaml, pocketResiduesDest)
        ## if specified in config, remove all other files and directories
        ## produced by fpocket/docking procedure
        if rmDirs:
            rmtree(runDir)
##########################################################################
##########################################################################

##########################################################################
def get_ligand_res_names(ligandDir: os.PathLike, ligands: list) -> list:
    """
    Simple loop function to look through ligand pdb files and extract three-letter residue names
    """
    resNames: list = []
    for ligand in ligands:
        ligandPdb: Union[PathLike, str] = p.join(ligandDir, f"{ligand}.pdb")
        ligandDf: pd.DataFrame = pdbUtils.pdb2df(ligandPdb)
        resName: str = ligandDf["RES_NAME"].tolist()[0]
        resNames.append(resName)

    return resNames
    









