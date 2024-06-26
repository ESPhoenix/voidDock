# import basic libraries
import os
import os.path as p
from subprocess import run, PIPE
from shutil import copy
# clean code
from typing import Union, Tuple, List
from os import PathLike



#######################################################################
def run_gnina(outDir: Union[PathLike, str],
              gninaConfig: Union[PathLike, str],
              gninaExe: Union[PathLike, str]) -> None:
      
    logFile: Union[PathLike, str] = p.join(outDir, "vina_docking.log")
    with open(logFile, "a") as logFile:
        run(f"{gninaExe} --config {gninaConfig}",
            shell=True, stdout=logFile, stderr=PIPE) ## TODO: make this work without shell=True
#######################################################################
def generate_gnina_flexible_residues():
    ...
    ##TODO: come back to flex docking when void portion is done


#######################################################################
def write_gnina_config(
        outDir: Union[PathLike, str],
        receptorPdb: Union[PathLike, str],
        ligands: List[Union[PathLike, str]],
        boxCenter: list,
        boxSize: list,
        flexibleResidueSyntax: str = None,
        exhaustiveness: int=16,
        numModes: int=10,
        cpus: int=2,
        seed: int=42,
        flex: bool=False) -> Tuple[Union[PathLike, str],  Union[PathLike, str]]:
    
  gninaConfigFile: Union[PathLike, str] = p.join(outDir, f"gnina_conf.txt")

  with open(gninaConfigFile, "w") as outFile:
    ## input pdbqt files, use flexible residues if required
    outFile.write(f"receptor = {receptorPdb}\n")
    if flex:  
        outFile.write(f"flexres = {flexibleResidueSyntax}\n\n")
    for ligandPdbqt in ligands:
        outFile.write(f"ligand = {ligandPdbqt}\n")

    ## docking box center coords (calculated as center of pocket)
    outFile.write(f"center_x = {str(boxCenter[0])}\n")
    outFile.write(f"center_y = {str(boxCenter[1])}\n")
    outFile.write(f"center_z = {str(boxCenter[2])}\n\n")
    ## docking box size
    outFile.write(f"size_x = {str(boxSize)}\n")
    outFile.write(f"size_y = {str(boxSize)}\n")
    outFile.write(f"size_z = {str(boxSize)}\n\n")
    ## exhaustiveness
    outFile.write(f"exhaustiveness = {str(exhaustiveness)}\n")
    ## num modes
    outFile.write(f"num_modes = {str(numModes)}\n\n")
    ## seed
    outFile.write(f"seed = {str(seed)}\n\n")
    ## cpus
    outFile.write(f"cpu = {str(cpus)}\n\n")
    ## output pdb file
    dockedPdb = p.join(outDir, "binding_poses.pdb")
    outFile.write(f"out = {dockedPdb}\n\n")


  return gninaConfigFile, dockedPdb
#######################################################################

def set_up_directory(outDir: Union[PathLike, str],
                      pathInfo: dict,
                        dockingOrder: dict) -> Tuple[str, Union[PathLike, str], list, Union[PathLike, str]]:
    # read protein pdb file, get name and make new dir for docking, copy over
    # protein pdb
    protName: str = dockingOrder["protein"]
    ligands: list = dockingOrder["ligands"]
    ligandDir: Union[PathLike, str] = pathInfo["ligandDir"]
    protDir: Union[PathLike, str] = pathInfo["protDir"]

    protPdb: Union[PathLike, str] = p.join(protDir, f"{protName}.pdb")

    # read ligand pdb and copy to new run directory
    ligandNames: list = []
    ligPdbqts:  list = []
    for ligandName in ligands:
        ligandNames.append(ligandName)
        ligPdbqt: Union[PathLike, str] = p.join(ligandDir, f"{ligandName}.pdbqt")
        ligPdbqts.append(ligPdbqt)   
    ## make a unique name for this docking simulation, 
    ## create directory with this name within outDir
    ligandTag: str = "_".join(ligandNames)
    runDir: Union[PathLike, str] = p.join(outDir, f"{protName}_{ligandTag}")
    os.makedirs(runDir, exist_ok=True)
    # copy over protPdb, move location of var
    protPdbDest: Union[PathLike, str] = p.join(runDir, f"{protName}.pdb")
    copy(protPdb, protPdbDest)
    protPdb: Union[PathLike, str] = protPdbDest

    return protName, protPdb, ligPdbqts, runDir
