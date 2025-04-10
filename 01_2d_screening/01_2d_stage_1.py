#!/usr/bin/env python
# coding: utf-8

'''
Written by: Sergei Kotelnikov
Modified by: Omeir Khan

Run a 2D screen on the Enamine database .cxsmiles
files using a given list of smarts patterns, provided
by --smarts_list
'''

# In[1]:


import os
import tqdm
import glob
import argparse
from rdkit import Chem
import multiprocessing as mp
from itertools import product
cpu_count = mp.cpu_count() - 1

# Entire "ENAMINE REAL" database
SMILES_FILES = ['/projectnb/docking/omeir/Enamine_REAL_db/Enamine_REAL_HAC_6_21_420M_CXSMILES.cxsmiles',
                '/projectnb/docking/omeir/Enamine_REAL_db/Enamine_REAL_HAC_22_23_471M_CXSMILES.cxsmiles',
                '/projectnb/docking/omeir/Enamine_REAL_db/Enamine_REAL_HAC_24_394M_CXSMILES.cxsmiles',
                '/projectnb/docking/omeir/Enamine_REAL_db/Enamine_REAL_HAC_25_557M_CXSMILES.cxsmiles',
                '/projectnb/docking/omeir/Enamine_REAL_db/Enamine_REAL_HAC_26_833M_Part_1_CXSMILES.cxsmiles',
                '/projectnb/docking/omeir/Enamine_REAL_db/Enamine_REAL_HAC_26_833M_Part_2_CXSMILES.cxsmiles',
                '/projectnb/docking/omeir/Enamine_REAL_db/Enamine_REAL_HAC_27_1.1B_Part_1_CXSMILES.cxsmiles', 
                '/projectnb/docking/omeir/Enamine_REAL_db/Enamine_REAL_HAC_27_1.1B_Part_2_CXSMILES.cxsmiles',
                '/projectnb/docking/omeir/Enamine_REAL_db/Enamine_REAL_HAC_28_1.2B_Part_1_CXSMILES.cxsmiles',
                '/projectnb/docking/omeir/Enamine_REAL_db/Enamine_REAL_HAC_28_1.2B_Part_2_CXSMILES.cxsmiles',
                '/projectnb/docking/omeir/Enamine_REAL_db/Enamine_REAL_HAC_29_38_988M_Part_1_CXSMILES.cxsmiles',
                '/projectnb/docking/omeir/Enamine_REAL_db/Enamine_REAL_HAC_29_38_988M_Part_2_CXSMILES.cxsmiles'
                ]

parser = argparse.ArgumentParser()
parser.add_argument('--smarts_list', '-l', help='A file containing smarts patterns to use in the 2D screen. Only one SMARTS pattern must be writte on each line. Substructure hits from the smarts pattern on the first line will be saved to a folder called {out_dir}/anchor_1/2d_stage.txt')
parser.add_argument('--out_dir', '-o', help='Name of the directory where outputs will be stored (default = ./results_2d_screen)', default='results_2d_screen/')

args = parser.parse_args()

def mp_func(mp_plan):
    smiles_file = mp_plan[0]
    smiles_file_name = os.path.basename(smiles_file)
    plan = mp_plan[1]
    output = {}
    for a in plan:
        output[a] = []

    with open(smiles_file) as f:
        f.readline()
        for line in f:
            smiles, idnumber = line.split('\t')[:2] # Assume SMILES and IDNumber are the first 2 cols
            mol = Chem.MolFromSmiles(smiles)
            
            # Save if there is a substructure match
            if mol:
                for name, patt_dict in plan.items():
                    smarts = patt_dict['smarts']
                    anchor_matches = mol.GetSubstructMatches(patt_dict['anchor_patt'])
                    tail_matches = []
                    
                    if len(anchor_matches) > 0:
                        output[name].append((smiles, idnumber, smarts))
                        print(smiles, idnumber, name, smiles_file_name)
                    
                    # Debugging
                    if len(output[name]) > 10:
                        return output 

    return output

def main():
    os.makedirs(args.out_dir, exist_ok=True)

    plan = {}
    with open(args.smarts_list) as f:
        idx = 1
        for line in f:
            smarts_patt = line.strip()

            patt = Chem.MolFromSmarts(smarts_patt)
            if patt:
                plan[f'anchor_{idx}'] = {'anchor_patt' : patt, 'smarts': smarts_patt}
                idx += 1
            else:
                print(f'RDKit cannot read {smarts_patt}!')
                

    print(f'Query plan: {plan}')
    
    # Parallel search of all database .cxsmiles files
    mp_plans = []
    for f in SMILES_FILES:
        mp_plans.append((f, plan))

    with mp.Pool(cpu_count) as pool:
        outputs = pool.map(mp_func, mp_plans)

    # Save results for each file
    for name in plan:
        if os.path.exists(f'{args.out_dir}/{name}') == False:
          os.mkdir(f'{args.out_dir}/{name}')

        file_lines = []
        for output in outputs:
            for smiles, idnumber, smarts in output[name]:
                file_lines.append(f'{smiles}\t{idnumber}\t{smarts}')
        #file_name = os.path.join(name, '2d_stage_1.txt')
        with open(f'{args.out_dir}/{name}/2d_stage_1.txt', 'w') as f:
            f.write('\n'.join(file_lines))
        print(len(file_lines))
        print(len(set(file_lines)))

if __name__=='__main__':
    main()
