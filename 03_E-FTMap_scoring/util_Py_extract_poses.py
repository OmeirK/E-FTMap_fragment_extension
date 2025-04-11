import argparse
import pandas as pd
from rdkit import Chem

parser =  argparse.ArgumentParser()

parser.add_argument('--eftmap_score_file', '-es', help='The .tsv file with E-FTMap scores (and paths to poses!)')
parser.add_argument('--n_ligands', '-n', help='Extract poses for the top N scoring ligands (default = 2000)', type=int, default=2000)
parser.add_argument('--outfile', '-o', help='The name of the sdf file with extractred poses (default = extracted_poses.sdf)', default='extracted_poses.sdf')

args = parser.parse_args()

def find_poses(eftmap_score_f, n_poses):
    df = pd.read_csv(eftmap_score_f, delimiter='\t')

    # In case the user gives an unsorted score file
    df_sort = df.sort_values(by=['Total'], ascending=False)

    data_dict = {}
    score_data = []
    for i, ligid in enumerate(df_sort['Lig_Name']):
        mol_idx = int(df_sort.iloc[i]['Mol_Idx'])
        sdf_f = df_sort.iloc[i]['Origin_File']

        # Debug:
        score = df_sort.iloc[i]['Total']
        
        data_dict[ligid] = {'sdf' : sdf_f,
                            'mol_idx': mol_idx}
        score_data.append((ligid, score))
        
    return data_dict, score_data

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
    print('[1/2] Reading scores...')
    pose_dict, score_data = find_poses(args.eftmap_score_file, args.n_ligands)

    print(f'[2/2] Extracting poses to {args.outfile}...')
    extract_poses(score_data, pose_dict, args.n_ligands, args.outfile)
       

if __name__=='__main__':
    main()
