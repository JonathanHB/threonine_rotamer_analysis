#Jonathan Borowsky
#Grabe lab
#2/4/26

#rotamer, hydrogen bonds, and secondary structure calculation for collaboration with Rutter lab
#############################################################################

#-------------------------------- imports -------------------------------

import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from MDAnalysis.analysis.dssp import DSSP

import compute_observable_on_trjs_x01_v3


def compute_rhs(trj_path, top_path, whichprotein, savefilename):

    upperpath = "/media/X01Raid01/Data_Backup/home/jborowsky/aac1-ucp1-aa-and-threonine/long-aac1-ucp1-processing/resi_data/"
    _savefilename = upperpath + savefilename

    #-------------------------------- load trajectory -------------------------------

    whichresidue = "ser"

    u = mda.Universe(top_path, trj_path, topology_format='ITP')

    n_frames = len(u.trajectory)
    n_residues = len(u.select_atoms(f"resname {whichresidue.upper()}").residues)


    #######################################################################################
    #                                 SECONDARY STRUCTURE

    #-------------------------------- calculate secondary structure for each frame -------------------------------

    secstr = DSSP(u).run().results.dssp
    secstr_num = np.where(secstr=='H', 1, 0)

    #-------------------------------- find threonines on helices -------------------------------
    #0-indexed to match secondary structure
    thr_indices = query_resis.resids-u.residues[0].resid

    secstr_thr = secstr_num[thr_indices]

    # helix_counts = np.zeros((n_residues, n_frames))
    # helix_assignments = np.zeros((n_residues, n_frames)).astype(int)

    # for i, th in enumerate(thr_indices):
    #     #the first few residues don't have enough residues on the n side of them to form another helix turn, and would throw index out of bounds errors
    #     if th>=4: 
    #         #select each threonine and the four preceding residues, which include the carbonyl of interest
    #         support_resis = secstr_num[:,th-4:th] 
    #         support_helical_count = np.sum(support_resis, axis = 1)
    #         helix_assignments[i] = np.where(support_helical_count==4, 1, 0)
    #         helix_counts[i] = support_helical_count


    # np.save(_savefilename+f"-{whichresidue}-secstr" + ".npy", helix_assignments.astype(np.uint8))

    np.save(_savefilename+f"-{whichresidue}-secstr_num" + ".npy", secstr_thr.astype(np.uint8))



localtest=False

if localtest:
    whichprotein = "ucp1"
    savefilename = "test"

    top_path = f"/home/jonathan/Documents/grabelab/rutter-lab-collaboration/test_trj/{whichprotein}/input/topol.top"
    trj_path = f"/home/jonathan/Documents/grabelab/rutter-lab-collaboration/test_trj/{whichprotein}/traj_comp.part0001-centered.xtc"
    compute_rhs(trj_path, top_path, whichprotein, savefilename)

else:
    compute_observable_on_trjs_x01_v3.compute_observable_x01(compute_rhs)
