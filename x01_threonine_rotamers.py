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

    whichresidue = "val"

    #note that these are the residue numbers, which are one larger than the corresponding indices
    threonines_ucp1 = [[193, 189, 121, 30, 122, 35],[226, 103, 97, 237, 259, 165, 166, 67, 62],[176, 177, 170, 156, 133, 60, 264],[11, 48, 108, 154]]
    threonines_aac1 = [[23,125,161,196,218,232],[83,138,220],[253],[154,251]]
    serines_ucp1 = [[18, 86, 89, 217, 279],[87, 183, 229],[114, 230, 273],[49, 50, 75, 111, 144, 242, 249, 296, 298]]
    serines_aac1 = [[],[21, 118, 126, 166, 227],[6,68,179,212,275],[1,41,45,241]]
    valines_ucp1 = [[282, 286, 216, 192, 184, 125, 127, 98, 28, 136, 235, 257, 64, 65],[38, 42, 138, 135, 202, 234],[211, 14, 271, 57, 232, 250],[51]]
    valines_aac1 = [[16, 31, 37, 57, 77, 180, 193, 219, 222, 226, 233, 277, 286],[25, 58, 130, 288],[178, 211, 254],[102, 144, 207, 297]]

    if whichresidue == "thr":
        if whichprotein == "aac1":
            residues_byloc = threonines_aac1
            offset = 1
        elif whichprotein == "ucp1":
            residues_byloc = threonines_ucp1
            offset = 9
        else:
            print(f"invalid protein {whichprotein}")
            return 0
        
    elif whichresidue == "ser":
        if whichprotein == "aac1":
            residues_byloc = serines_aac1
            offset = 1
        elif whichprotein == "ucp1":
            residues_byloc = serines_ucp1
            offset = 9
        else:
            print(f"invalid protein {whichprotein}")
            return 0
    
    elif whichresidue == "val":
        if whichprotein == "aac1":
            residues_byloc = valines_aac1
            offset = 1
        elif whichprotein == "ucp1":
            residues_byloc = valines_ucp1
            offset = 9
        else:
            print(f"invalid protein {whichprotein}")
            return 0

    else:
        print(f"invalid protein {whichresidue}")
        return 0   

    
    #residue indices matching mdanalysis object (0-indexed rather than using the PDB numbering)
    residue_inds_byloc = [[ii-offset for ii in i] for i in residues_byloc]

    #determine array dimensions necessary to encode residues in all location types in a single rectilinear numpy array
    # this array is padded with values of 400 (outside the angle range returned by the dihedral calculation) 
    # where there are not enough residues in a particular location type
    n_frames = len(u.trajectory)
    n_loctypes = len(residues_byloc)
    n_resis_per_loctype = max(len(lst) for lst in residues_byloc)


    #-------------------------------- verify residue list consistency -------------------------------

    flat_resis = np.array([si for s in residue_inds_byloc for si in s])
    flat_resis.flatten()

    #check for duplicated residue indices in the list
    duplicates = (len(flat_resis) != len(np.unique(flat_resis)))
    if duplicates:
        print("error: residue list includes residues multiple times")
        return 0

    #check for residues indices in the list corresponding to the wrong residue in the protein structure
    all_correct_resn = np.all([res.resname == whichresidue.upper() for res in protein_trj.residues[flat_resis]])
    if not all_correct_resn:
        print(f"error: residue list contains residues other than {whichresidue}")
        return 0

    #check for residues present in the protein structure but not the list
    total_resis = len([res.resname for res in protein_trj.residues if res.resname.lower() == whichresidue])
    if total_resis > len(np.unique(flat_resis)):
        print(f"error: not all {whichresidue} residues are listed")
        return 0
    elif total_resis < len(np.unique(flat_resis)):
        print(f"error: more residues are listed than exist in the protein")
        return 0


    #-------------------------------- calculate dihedrals -------------------------------

    #array to hold results for all threonines 
    #of shape [number of location types, maximum number of residues at any one location type, number of frames]
    angles_all = np.ones((n_loctypes, n_resis_per_loctype, n_frames))*400

    for loc_ind in range(n_loctypes):
        resi_inds = residue_inds_byloc[loc_ind]
        if resi_inds == []:
            continue
        #see https://userguide.mdanalysis.org/1.1.1/examples/analysis/structure/dihedrals.html
        chi1_atoms = [res.chi1_selection() for res in protein_trj.residues[resi_inds]]
        chi1_dihedrals = dihedrals.Dihedral(chi1_atoms).run()
        dihedral_angles = chi1_dihedrals.results.angles
        angles_all[loc_ind, 0:len(resi_inds)] = dihedral_angles.transpose()


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

    np.save(_savefilename+f"-{whichresidue}-rotamer-trj-discrete" + ".npy", rotamer_trj.astype(np.uint8))


    #-------------------------------- count rotamter occurrences -------------------------------
    # note that due to trajectory length variability frequencies are not to be calculated until later
    rotamer_counts = np.zeros((n_loctypes, n_resis_per_loctype, n_rotamers))

    for rotamer_ind in range(n_rotamers):
        rotamers_onehot = np.where(rotamer_trj == rotamer_ind, 1, 0)
        rotamer_counts[:,:,rotamer_ind] = np.sum(rotamers_onehot, axis = 2)

    np.save(_savefilename+f"-{whichresidue}-rotamer-counts" + ".npy", rotamer_counts)


compute_observable_on_trjs_x01_v3.compute_observable_x01(compute_rotamers)
