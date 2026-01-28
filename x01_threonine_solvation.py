import numpy as np
import itertools

import matplotlib.pyplot as plt

import mdtraj as md

import time
import os
import sys

sys.path.insert(1, f'{os.getcwd()}/../utility')
import compute_observable_on_trjs_x01_v3


def get_trj_contacts(trj, protein_query, ligand_query, contact_threshold):
    
    #select protein and ligand atoms
    indices1 = trj.top.select(ligand_query)
    indices2 = trj.top.select(protein_query)

    # list of atom pairs
    index_pairs = list(itertools.product(indices1, indices2))
    
    # calculate pairwise distances
    # periodic = True is an issue for rcsb pdb structures with unit cell information if the residues are far enough apart
    distances = md.compute_distances(trj, index_pairs, periodic=True, opt=True).reshape(trj.n_frames, len(indices1), len(indices2))

    #calculate the minimum distance from each protein atom to any ligand atom
    #mindists = np.min(distances, axis=1)
    contacts = np.where(distances < contact_threshold, 1, 0)
    print(contacts.shape)
    n_contacts = np.sum(contacts, axis=1)

    return n_contacts


def thr_solv(trj_path, top_path, protein, savefilename):

    #TODO comment after processing AAC1 data
    if protein != "aac1":
        return 0

    upperpath = "/media/X01Raid01/Data_Backup/home/jborowsky/aac1-ucp1-aa-and-threonine/long-aac1-ucp1-processing/thr_solvation/"
    _savefilename = upperpath + savefilename

    print(f"trajectory = {trj_path}")
    print(f"topology = {top_path}")

    #parameters
    contact_threshold = 0.33 #nm
    water_query = "resname HOH and element O"
    
    #residue categories
    #on a helix facing the membrane core     <--,----- these two have helix carbonyls to hydrogen bond with; the two other categories do not
    #on a helix facing water or headgroups   <--'
    #on a helix N terminal (variable solvation)
    #on a loop (all of these are in water)
    threonines_bytype_ucp1 = [[193, 189, 121, 30, 122, 35],[226, 103, 97, 237, 259, 165, 166, 67, 62],[176, 177, 170, 156, 133, 60, 264],[11, 48, 108, 154]]
    threonines_bytype_aac1 = [[23,125,161,196,218,232],[83,138,220],[253],[154,251]]

    if protein == "aac1":
        threonines_bytype = threonines_bytype_aac1
    elif protein == "ucp1":
        threonines_bytype = threonines_bytype_ucp1
    else:
        print(f"invalid protein {protein}")
        return 0

    #queries
    threonines_flat = [i for sublist in threonines_bytype for i in sublist]
    thrresi = " or ".join(["resSeq "+str(i) for i in threonines_flat])
    thr_query = f"protein and ({thrresi}) and name OG1"

    #load trajectory
    ref_trj = md.load(top_path)
    ref_inds = ref_trj.top.select("not element H")
    trj = md.load(trj_path, top = top_path, atom_indices = ref_inds)

    print(f"trajectory has {trj.n_frames} frames")

    #calculate solvation
    n_contacts = get_trj_contacts(trj[::10], thr_query, water_query, contact_threshold)

    nsolv = n_contacts.astype(np.uint8)

    np.save(_savefilename+"-nsolv-s10" + ".npy", nsolv)



compute_observable_on_trjs_x01_v3.compute_observable_x01(thr_solv)