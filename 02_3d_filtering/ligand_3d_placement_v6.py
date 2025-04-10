'''
Written by: Omeir Khan/Sergei Kotelnikov

This script will generate 3D conformers for a given molecule,
align it to a substructure specified by the --anchor_smarts,
and then remove conformers that clash with the receptor atoms
in --rec_pdb

v3: [X] Cluster poses and only save cluster centers
v4: [X] Uhhhh I forgot to write what i changed :S
v5: [X] Add an RMSD filter for poses that were generated
    [X] Add a parameter for number of conformers to generate
v6: [] Disable multi-receptor docking to simplify the inputs
       for public use
'''


import os
import glob
import argparse
import contextlib
import numpy as np
import prody as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.ML.Cluster import Butina

parser = argparse.ArgumentParser()
parser.add_argument('--rec_pdb', '-r', help='Path to the receptor PDB file you want to dock to. Ensure that it is in the same frame of reference as the --template_lig that is provided')
parser.add_argument('--template_lig', '-tl', help='Path to a .mol file containing coordinates of the template ligand/fragment to use for docking')
parser.add_argument('--anchor_smarts', '-as', help='SMARTS pattern for the anchor')
parser.add_argument('--lig_smiles', '-ls', help='SMILES string for the ligand you want to dock')
parser.add_argument('--lig_name', '-ln', help='Name of the ligand you are docking (default=lig)', default='lig')
parser.add_argument('--outdir', '-od', help='Path to a directory to store the output in')
parser.add_argument('--noRandCoord', '-nrc', help='Enable this flag to disable use of the "randCoord" parameter for EKTDG conformer generation. This will place substructure atoms in the EXACT position of the template ligand, but at a significant computational cost.', action='store_false', default=True)
parser.add_argument('--n_confs', '-nc', help='The number of conformers to generate for template-based docking (default=100)', default=100, type=int)
parser.add_argument('--rmsd_cutoff', '-rc', help='The maximum RMSD a docked pose must have to the template ligand (default = 1.0)', default=1.0, type=float)
parser.add_argument('--clash_cutoff', '-cc', help='The minimum distance a docked pose can be from the receptor in order to not clash (default = 1.5)', default=1.5, type=float)

args = parser.parse_args()

@contextlib.contextmanager
def cwd(path):
    oldpwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(oldpwd)

# Map atoms in the target ligand to atoms in the template
def cs_sym_mappings(target_mol, template_mol, cs_smarts): # accounts for target_mol symmetry
    cs_patt = Chem.MolFromSmarts(cs_smarts)
    target_cs_matches = target_mol.GetSubstructMatches(cs_patt, uniquify=False)
    template_cs_matches = template_mol.GetSubstructMatches(cs_patt, uniquify=False)

    # Debugging #
    #print(target_cs_matches)
    #print(template_cs_matches)
    #print(Chem.MolToSmiles(target_mol), Chem.MolToSmiles(template_mol), cs_smarts)
    #print(template_mol.HasSubstructMatch(cs_patt))
    
    mappings = set()
    for target_cs_match in target_cs_matches:
        for template_cs_match in template_cs_matches:
            mapping = tuple(sorted(zip(target_cs_match, template_cs_match), key=lambda x: x[1]))
            mappings.add(mapping)
    
    mol_sym_matches = target_mol.GetSubstructMatches(target_mol, uniquify=False)
    mappings_reduced = []
    while len(mappings) > 0:
        mapping = list(mappings.pop())
        mappings_reduced.append(mapping)
        redundant_mappings = {tuple((mol_sym_match[i], j) for i, j in mapping) for mol_sym_match in mol_sym_matches}
        mappings -= redundant_mappings
    return mappings_reduced

# Calculate minimum RMSD between two sets of coordinates
def calc_min_rmsd(coords1, coords2, mappings):
    rmsds = []
    for mapping in mappings:
        match1, match2 = map(list, zip(*mapping))
        diff = coords2[match2] - coords1[match1]
        rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        rmsds.append(rmsd)
    return min(rmsds)

# Cluster a given set of conformers for a molecule
def clusterization(mol, cid_list, cluster_radius):
    sym_mappings = [list(enumerate(match)) for match in mol.GetSubstructMatches(mol, uniquify=False)]
    rmsmat = []
    size = len(cid_list)
    for i in range(size):
        coordsi = mol.GetConformer(cid_list[i]).GetPositions()
        for j in range(i):
            coordsj = mol.GetConformer(cid_list[j]).GetPositions()
            rmsmat.append(calc_min_rmsd(coordsi, coordsj, sym_mappings))
    clusters = Butina.ClusterData(rmsmat, size, cluster_radius, isDistData=True, reordering=True)
    cluster_centers_cids = [cid_list[cluster[0]] for cluster in clusters]
    return cluster_centers_cids

rec_pdb_files = [args.rec_pdb]
ref_mol_files = [args.template_lig]
smarts = args.anchor_smarts
smiles = args.lig_smiles.split('.')[0]
clash_dist = args.clash_cutoff
rmsd_cutoff = args.rmsd_cutoff

mol_ = Chem.MolFromSmiles(smiles) # Smiles should be correctly protonated!
for rec_pdb_file in rec_pdb_files: # I guess we should keep minimizations with different receptors separately (unlike anchors: ref_mols and mappings), they give different energy baselines
    rec_ag = pd.parsePDB(rec_pdb_file) # without hydrogens (we will add hydrogens during minimization), without sequence gaps, missing tails are ok (unless they are near the binding pocket)
    rec_coords = rec_ag.select('protein').getCoords() # you can change the selection depending on what you want to be present in receptor for heavy atom clash check
    rec_pdb_name = os.path.basename(rec_pdb_file)
    rec_pdbid = rec_pdb_name[:4]

    for ref_mol_file in ref_mol_files:
        ref_mol = Chem.MolFromMolFile(ref_mol_file)
        ref_mol_name = os.path.basename(ref_mol_file)

        sym_mappings = cs_sym_mappings(mol_, ref_mol, smarts) # or you can write manually your favorite mapping
        print(f'Atom mapping from {ref_mol_name} to {args.lig_name} ({args.lig_smiles}) ...')
        print(f'\t{len(sym_mappings)} mapping(s) have been generated')
        #print(sym_mappings)
        for mi, mapping in enumerate(sym_mappings):
            cmap = {i: ref_mol.GetConformer().GetAtomPosition(j) for i, j in mapping}
            sdf_file = f'{args.outdir}/{args.lig_name}_{rec_pdb_name[:-4]}.{ref_mol_name[:-4]}.{mi}.sdf'

            if os.path.exists(sdf_file):
                print(f'\t{sdf_file} already exists! Skipping it.')
                continue                
            
            print(f'\tFinding poses for mapping {mi}...')
            #print(mapping)
            mol = Chem.AddHs(mol_)
            mol.SetProp('_Name', args.lig_name)
            cids = AllChem.EmbedMultipleConfs(mol, numConfs=args.n_confs, coordMap=cmap, numThreads=1, useRandomCoords=args.noRandCoord)
            
            print(f'\t\t{len(cids)} poses generated')
            # Skip if no conformers were generated
            if len(cids) == 0:
                continue

            print(f'\t\tChecking poses for clash with receptor {rec_pdb_name}...')
            noclash_cids = []
            for cid in cids:
                rmsd = rdMolAlign.AlignMol(mol, ref_mol, prbCid=cid, atomMap=mapping)
                if rmsd <= rmsd_cutoff:
                    coords = Chem.RemoveHs(mol).GetConformer(cid).GetPositions() # ligand heavy atoms
                    dist = pd.buildDistMatrix(coords, rec_coords).min()
                    if dist >= clash_dist:
                        #print(f'\t\t\tRMSD {rmsd:.2f}', rmsd_cutoff)
                        #print(f'\t\t\tDist {dist:.2f}', clash_dist)
                        noclash_cids.append(cid)
            
            print(f'\t\t\t{len(noclash_cids)} conformers pass RMSD and clash filters!')

            # Cluster remaining conformers with an RMSD cutoff of 0.5
            cluster_center_cids = clusterization(mol, noclash_cids, 0.5)
            
            print(f'\t\tSaving {len(cluster_center_cids)} clusters...')
            with Chem.SDWriter(sdf_file) as w:
                for cid in cluster_center_cids:
                    w.write(mol, confId=cid)


            # Check if anything was written to the sdf_file
            # Delete it if it is empty
            if os.path.exists(sdf_file):
                with open(sdf_file) as f:
                    sdf_lines = f.readlines()

                if len(sdf_lines) == 0:
                    os.remove(sdf_file)
