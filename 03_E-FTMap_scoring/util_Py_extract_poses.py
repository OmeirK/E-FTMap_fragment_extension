import argparse
import pandas as pd
from rdkit import Chem

parser =  argparse.ArgumentParser()

parser.add_argument('--combo_score_file', '-cs', help='The .tsv file with E-FTMap, Diffusion, and Combo scores')
parser.add_argument('--eftmap_score_file', '-es', help='The .tsv file with E-FTMap scores only (containing paths to poses)')
parser.add_argument('--n_ligands', '-n', help='Extract poses for the top N scoring ligands (default = 2000)', type=int, default=2000)
parser.add_argument('--outfile', '-o', help='The name of the sdf file with extractred poses (default = extracted_poses.sdf)', default='extracted_poses.sdf')
parser.add_argument('--eftmap', '-e', action='store_true', default=False)
parser.add_argument('--diffusion', '-d', action='store_true', default=False)
parser.add_argument('--combo', '-c', action='store_true', default=False)

args = parser.parse_args()

def read_combo_scores(combo_f, eftmap=False, diffusion=False, combo=False):
    df = pd.read_csv(combo_f, delimiter='\t')
    data_l = []
    for i, ligid in enumerate(df['Enamine_ID']):
        if eftmap == True:
            score = df['E-FTMap_Score'][i]
        if diffusion == True:
            score = df['Diff_Conf'][i]
        if combo == True:
            score = df['Projected_Score'][i]
        
        data_l.append((ligid, score))

    # Sort the list by score (larger is better for all three metrics)
    data_l.sort(key=lambda x: x[1], reverse=True)

    return data_l

def find_poses(eftmap_score_f):
    df = pd.read_csv(eftmap_score_f, delimiter='\t')

    data_dict = {}
    read_sdfs = {}
    for i, ligid in enumerate(df['Lig_Name']):
        mol_idx = int(df['Mol_Idx'][i])
        sdf_f = df['Origin_File'][i]
        
        data_dict[ligid] = {'sdf' : sdf_f,
                            'mol_idx': mol_idx}

    return data_dict

def extract_poses(score_data, pose_dict, n_poses, outfile):
    read_sdfs = {}

    with Chem.SDWriter(outfile) as w:
        for i in range(0, n_poses):
            ligid = score_data[i][0]
            score = score_data[i][1]
            sdf_f = pose_dict[ligid]['sdf']
            pose_idx = pose_dict[ligid]['mol_idx']

            if sdf_f not in read_sdfs:
                print(f'\tLoading mols from {sdf_f}...')
                mols = Chem.SDMolSupplier(sdf_f)
                read_sdfs[sdf_f] = mols
            
            mol = read_sdfs[sdf_f][pose_idx]
            mol.SetProp('score', f'{score:.2f}')

            w.write(mol)

    print(f'{n_poses} ligands successfully extracted to {outfile}!')


def main():
    if [args.eftmap, args.diffusion, args.combo].count(True) > 1:
        raise ValueError('Only enable one of the following: --eftmap, --diffusion, --combo')
    print([args.eftmap, args.diffusion, args.combo])

    if args.eftmap:
        print('Reading eftmap scores!')
    if args.diffusion:
        print('Reading diffusion scores!')
    if args.combo:
        print('Reading projected combo scores!')



    print('[1/3] Reading scores...')
    score_data = read_combo_scores(args.combo_score_file, eftmap=args.eftmap, diffusion=args.diffusion, combo=args.combo)
    print(score_data[:10])

    print('[2/3] Finding pose files...')
    pose_dict = find_poses(args.eftmap_score_file)
    print(f'[3/3] Extracting poses to {args.outfile}...')
    extract_poses(score_data, pose_dict, args.n_ligands, args.outfile)
       

if __name__=='__main__':
    main()
