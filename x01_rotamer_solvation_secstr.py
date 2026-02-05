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

    upperpath = "/media/X01Raid01/Data_Backup/home/jborowsky/aac1-ucp1-aa-and-threonine/long-aac1-ucp1-processing/thr_rotamers_v2/"
    _savefilename = upperpath + savefilename

    #-------------------------------- load trajectory -------------------------------

    whichresidue = "thr"

    u = mda.Universe(top_path, trj_path, topology_format='ITP')

    n_frames = len(u.trajectory)
    n_residues = len(u.select_atoms(f"resname {whichresidue.upper()}").residues)


    #######################################################################################
    #                                 HYDROGEN BONDS

    #-------------------------------- calculate hydrogen bonds -------------------------------

    donornames = {"thr":"HG1", "ser":"HG1"}
    acceptor_names = {"thr":"OG1", "ser":"OG"}

    d_a_cutoff=3.3

    h1 = HBA(universe=u, hydrogens_sel=f"resname {whichresidue.upper()} and name {donornames[whichresidue]}",
            acceptors_sel="resname TIP3 and name OH2", update_selections=False, d_a_cutoff=d_a_cutoff)
    h1.run()

    h2 = HBA(universe=u, hydrogens_sel="resname TIP3 and name H1 H2", 
            acceptors_sel=f"resname {whichresidue.upper()} and name {acceptor_names[whichresidue]}", update_selections=False, d_a_cutoff=d_a_cutoff)
    h2.run()


    #-------------------------------- count hydrogen bonds by residue -------------------------------

    #number of hydrogen bonds by frame and residue
    hbonds_all = np.zeros((n_residues, n_frames)).astype(int)

    query_resis = u.select_atoms(f"resname {whichresidue.upper()}").residues
    #donor_inds = [[a.index for a in r.atoms if a.name == "HG1"][0] for r in query_resis]
    acceptor_inds = [[a.index for a in r.atoms if a.name == acceptor_names[whichresidue]][0] for r in query_resis] #these are atom indices

    for i, ai in enumerate(acceptor_inds):

        times_1 = h1.results.hbonds[np.where(h1.results.hbonds[:,1].astype(int)==ai),0].astype(int)
        for t1 in times_1:
            hbonds_all[i, t1] += 1 

        times_2 = h2.results.hbonds[np.where(h2.results.hbonds[:,3].astype(int)==ai),0].astype(int)
        for t2 in times_2:
            hbonds_all[i, t2] += 1 


    np.save(_savefilename+f"-{whichresidue}-hbonds" + ".npy", hbonds_all.astype(np.uint8))


    #######################################################################################
    #                                 ROTAMERS

    #-------------------------------- calculate dihedrals -------------------------------

    #see https://userguide.mdanalysis.org/1.1.1/examples/analysis/structure/dihedrals.html
    chi1_atoms = [res.chi1_selection() for res in query_resis]
    chi1_dihedrals = dihedrals.Dihedral(chi1_atoms).run()
    dihedral_angles = chi1_dihedrals.results.angles
    x1_angles = dihedral_angles.transpose()

    #-------------------------------- discretize trajectory -------------------------------
    #mapping: 
    #         a < -120 : 0
    # -120 <= a <    0 : 1
    #    0 <= a <  120 : 2
    #  120 <= a         : 3

    rotamer_bounds = [-120, 0, 120]
    n_rotamers = len(rotamer_bounds)

    rotamer_trj = np.digitize(x1_angles, rotamer_bounds) #this splits the rotamer centered on +-180 degrees into two bins
    rotamer_trj = np.where(rotamer_trj == n_rotamers, 0, rotamer_trj).astype(int) #merge the bins for the rotamer at 180 degrees back together


    np.save(_savefilename+f"-{whichresidue}-rotamers" + ".npy", rotamer_trj.astype(np.uint8))


    #######################################################################################
    #                                 SECONDARY STRUCTURE

    #-------------------------------- calculate secondary structure for each frame -------------------------------

    secstr = DSSP(u).run().results.dssp
    secstr_num = np.where(secstr=='H', 1, 0)# + np.where(secstr=='E', 2, 0) #there are vanishingly few sheet assignments, which are all standalone misclassified loops


    #-------------------------------- find threonines on helices long enough to support backbone H bonding -------------------------------
    #0-indexed to match secondary structure
    thr_indices = query_resis.resids-u.residues[0].resid

    helix_counts = np.zeros((n_residues, n_frames))
    helix_assignments = np.zeros((n_residues, n_frames)).astype(int)

    for i, th in enumerate(thr_indices):
        #the first few residues don't have enough residues on the n side of them to form another helix turn, and would throw index out of bounds errors
        if th>=4: 
            #select each threonine and the four preceding residues, which include the carbonyl of interest
            support_resis = secstr_num[:,th-4:th] 
            support_helical_count = np.sum(support_resis, axis = 1)
            helix_assignments[i] = np.where(support_helical_count==4, 1, 0)
            helix_counts[i] = support_helical_count


    np.save(_savefilename+f"-{whichresidue}-secstr" + ".npy", helix_assignments.astype(np.uint8))



localtest=False

if localtest:
    whichprotein = "ucp1"
    savefilename = "test"

    top_path = f"/home/jonathan/Documents/grabelab/rutter-lab-collaboration/test_trj/{whichprotein}/input/topol.top"
    trj_path = f"/home/jonathan/Documents/grabelab/rutter-lab-collaboration/test_trj/{whichprotein}/traj_comp.part0001-centered.xtc"
    compute_rhs(trj_path, top_path, whichprotein, savefilename)

else:
    compute_observable_on_trjs_x01_v3.compute_observable_x01(compute_rhs)
