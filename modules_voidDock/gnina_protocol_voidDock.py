# import basic libraries
import os
import os.path as p
from subprocess import run
from shutil import copy
# inport util scripts
from modules_voidDock import shared_prep_voidDock as sharedPrep
from modules_voidDock import gnina_prep_voidDock as gninaPrep
# clean code
from typing import Union, Tuple
from os import PathLike

#######################################################################
def docking_protocol(config: dict, dockingOrder:dict) -> None:
    # read config
    pathInfo: Union[PathLike, str] = config["pathInfo"]
    outDir: Union[PathLike, str] = pathInfo["outDir"]
    cpusPerRun: Union[PathLike, str] = config["cpuInfo"]["cpusPerRun"]
    gninaExe: Union[PathLike, str] = config["methodInfo"]["gninaExe"]

    print()


    # set up run directory and output key variables
    protName, protPdb, ligPdbqts, runDir = gninaPrep.set_up_directory(outDir=outDir,
                                                             pathInfo=pathInfo,
                                                               dockingOrder=dockingOrder)

    # Use fpocket to identify correct pocket, calclate box center and return
    # residues in pocket
    targetPocketResidues: list = dockingOrder["pocketResidues"]
    boxCenter, pocketResidues = sharedPrep.directed_fpocket(protName=protName,
                                                 runDir=runDir,
                                                 pdbFile=protPdb,
                                                 targetPocketResidues=targetPocketResidues)


    if dockingOrder["mutatePocketToAla"]: 
        # Replace pocket residues with alanine
        protPdb: Union[PathLike, str] = sharedPrep.pocket_residues_to_alainine(protName=protName,
                                            pdbFile=protPdb,
                                              residuesToAlanine=pocketResidues,
                                                dockingOrder=dockingOrder,
                                                  outDir=runDir)


    flexibleDocking: bool = False
    if "flexibleResidues" in dockingOrder:
        if len(dockingOrder["flexibleResidues"]) > 0:
            flexibleResidueSyntax: str = gninaPrep.generate_gnina_flexible_residues()
            
            gninaConfig, dockedPdb = gninaPrep.write_gnina_config(outDir=runDir,
                                                receptorPdb=protPdb,
                                                ligands= ligPdbqts, 
                                                flexibleResidueSyntax = flexibleResidueSyntax,
                                                boxCenter=boxCenter,
                                                boxSize=30,
                                                cpus=cpusPerRun,
                                                flex=True)
    else:
        gninaConfig, dockedPdb = gninaPrep.write_gnina_config(outDir=runDir,
                                                receptorPdb=protPdb,
                                                ligands= ligPdbqts, 
                                                boxCenter=boxCenter,
                                                boxSize=30,
                                                cpus=cpusPerRun)
        

    gninaPrep.run_gnina(gninaConfig=gninaConfig,
              gninaExe=gninaExe,
              outDir=runDir)
