#Jonathan Borowsky
#Grabe lab
#01/28/26

#rotamer calculation for collaboration with Rutter lab
#############################################################################

#-------------------------------- imports -------------------------------

import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals

import sys
import os

#sys.path.insert(1, f'{os.getcwd()}/../utility')
import compute_observable_on_trjs_x01_v3


def compute_rotamers(trj_path, top_path, whichprotein, savefilename):

    #-------------------------------- load trajectory -------------------------------

    #toppath = "/home/jonathan/Documents/grabelab/rutter-lab-collaboration/test_trj/input/seg_0035.gro"
    #trjpath = "/home/jonathan/Documents/grabelab/rutter-lab-collaboration/test_trj/traj_comp.part0001-centered.xtc"
    u = mda.Universe(top_path, trj_path)
    protein_trj = u.select_atoms('protein')

    upperpath = "/media/X01Raid01/Data_Backup/home/jborowsky/aac1-ucp1-aa-and-threonine/long-aac1-ucp1-processing/thr_rotamers/"
    _savefilename = upperpath + savefilename


    #-------------------------------- select residues -------------------------------

    #whichprotein = "aac1"
    location_type_info = [
        "on a helix facing the membrane core",
        "on a helix facing water or headgroups",
        "on a helix N terminal with no helix carbonyl to hydrogen bond with (variable solvation)",
        "on a loop (all of these are in water)"
    ]

    #note that these are the residue numbers, which are one larger than the corresponding indices
    threonines_ucp1 = [[193, 189, 121, 30, 122, 35],[226, 103, 97, 237, 259, 165, 166, 67, 62],[176, 177, 170, 156, 133, 60, 264],[11, 48, 108, 154]]
    threonines_aac1 = [[23,125,161,196,218,232],[83,138,220],[253],[154,251]]

    if whichprotein == "aac1":
        threonines_byloc = threonines_aac1
        offset = 1
    elif whichprotein == "ucp1":
        threonines_byloc = threonines_ucp1
        offset = 9
    else:
        print(f"invalid protein {whichprotein}")
        #return 0

    threonine_inds_byloc = [[ii-offset for ii in i] for i in threonines_byloc]

    #determine array dimensions necessary to encode residues in all location types in a single rectilinear numpy array
    # this array is padded with values of 400 (outside the angle range returned by the dihedral calculation) 
    # where there are not enough residues in a particular location type
    n_frames = len(u.trajectory)
    n_loctypes = len(threonines_byloc)
    n_resis_per_loctype = max(len(lst) for lst in threonines_byloc)


    #-------------------------------- calculate dihedrals -------------------------------

    #array to hold results for all threonines 
    #of shape [number of location types, maximum number of residues at any one location type, number of frames]
    angles_all = np.ones((n_loctypes, n_resis_per_loctype, n_frames))*400

    for loc_ind in range(n_loctypes):
        thr_inds = threonine_inds_byloc[loc_ind]
        #see https://userguide.mdanalysis.org/1.1.1/examples/analysis/structure/dihedrals.html
        chi1_atoms = [res.chi1_selection() for res in protein_trj.residues[thr_inds]]
        chi1_dihedrals = dihedrals.Dihedral(chi1_atoms).run()
        dihedral_angles = chi1_dihedrals.results.angles
        angles_all[loc_ind, 0:len(thr_inds)] = dihedral_angles.transpose()


    #-------------------------------- discretize trajectory -------------------------------
    #mapping: 
    #        a <= -120 : 0
    # -120 < a <=    0 : 1
    #    0 < a <=  120 : 2
    #  120 < a <=  180 : 3
    #  180 < a         : 4 #used to represent empty array entries (i.e. at location types with fewer residues)
    rotamer_bounds = [-120, 0, 120, 180]
    n_rotamers = len(rotamer_bounds)-1

    rotamer_trj = np.digitize(angles_all, rotamer_bounds, right=True) #this splits the rotamer centered on +-180 degrees into two bins
    rotamer_trj = np.where(rotamer_trj == n_rotamers, 0, rotamer_trj) #merge the bins for the rotamer at 180 degrees back together

    np.save(_savefilename+"-thr-rotamer-trj-discrete" + ".npy", rotamer_trj.astype(np.uint8))


    #-------------------------------- count rotamter occurrences -------------------------------
    # note that due to trajectory length variability frequencies are not to be calculated until later
    rotamer_counts = np.zeros((n_loctypes, n_resis_per_loctype, n_rotamers))

    for rotamer_ind in range(n_rotamers):
        rotamers_onehot = np.where(rotamer_trj == rotamer_ind, 1, 0)
        rotamer_counts[:,:,rotamer_ind] = np.sum(rotamers_onehot, axis = 2)

    np.save(_savefilename+"-thr-rotamer-counts" + ".npy", rotamer_counts)


compute_observable_on_trjs_x01_v3.compute_observable_x01(compute_rotamers)
