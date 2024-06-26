

# clean code
from typing import Union
from os import PathLike
## voidDock libraries
from modules_voidDock import shared_prep_voidDock as sharedPrep
from modules_voidDock import vina_prep_voidDock as vinaPrep


def docking_protocol(config: dict, dockingOrder:dict) -> None:
    # read config
    pathInfo: Union[PathLike, str] = config["pathInfo"]
    outDir: Union[PathLike, str] = pathInfo["outDir"]
    cpusPerRun: Union[PathLike, str] = config["cpuInfo"]["cpusPerRun"]


    # set up run directory and output key variables
    protName, protPdb, ligPdbqts, runDir = vinaPrep.set_up_directory(outDir=outDir,
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
            flexibleDocking  = True
            rigidPdbqt, flexPdbqt = vinaPrep.gen_flex_pdbqts(protPdb = protPdb,
                                                    flexibleResidues = dockingOrder["flexibleResidues"],
                                                    outDir = runDir)

            # Write a config file for vina
            vinaConfig, dockedPdbqt = vinaPrep.write_vina_config(outDir=runDir,
                                                receptorPdbqt=rigidPdbqt,
                                                flexPdbqt = flexPdbqt,
                                                boxCenter=boxCenter,
                                                boxSize=30,
                                                cpus=cpusPerRun,
                                                flex=True)
                                    
    if not flexibleDocking:
        # Convert alanine PDB to PDBQT
        rigidPdbqt: Union[PathLike, str] = vinaPrep.pdb_to_pdbqt(inPdb=protPdb,
                                outDir=runDir,
                                jobType="rigid")

        # Write a config file for vina
        vinaConfig, dockedPdbqt = vinaPrep.write_vina_config(outDir=runDir,
                                            receptorPdbqt=rigidPdbqt,
                                            boxCenter=boxCenter,
                                            boxSize=30,
                                            cpus=cpusPerRun,
                                            flex=False)


    # Run vina docking
    vinaPrep.run_vina(outDir=runDir,
             configFile=vinaConfig,
             ligPdbqts=ligPdbqts)

        
    vinaPrep.process_vina_results(outDir= runDir,
                        dockedPdbqt = dockedPdbqt,
                        receptorPdbqt = rigidPdbqt,
                        dockingOrder = dockingOrder,
                        ligandDir= config["pathInfo"]["ligandDir"]
                        )
