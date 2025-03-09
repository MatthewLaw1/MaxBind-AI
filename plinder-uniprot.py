import pandas as pd
import ast
import json
from rdkit import Chem
from rdkit.Chem import AllChem

# Read the first CSV file (adjust the file path as needed)
df1 = pd.read_csv('/Users/ak87/ownCloud/owncloud/Calculations/HACKATHON/plindex_single_proper.csv')

# Read the second file. If it's tab-delimited, specify the separator; 
# otherwise, adjust accordingly.
df2 = pd.read_csv('/Users/ak87/ownCloud/owncloud/Calculations/HACKATHON/pdb_chain_uniprot.tsv', sep='\t', skiprows=1)

# Define a helper function to extract the chain letter from the list string.
def extract_chain(chain_str):
    """
    Given a string like "['1.A']", this function extracts the chain part (i.e., "A").
    It safely evaluates the string as a list and splits the element by '.'.
    """
    try:
        # Safely convert the string representation of a list to an actual list
        chains = ast.literal_eval(chain_str)
        # Assume the list has one element, e.g., "1.A"
        chain_value = chains[0]
        # Split on the dot to extract the chain identifier
        if '.' in chain_value:
            return chain_value.split('.')[1]
        else:
            return chain_value
    except Exception as e:
        # In case of error, return None (or you can choose to handle it differently)
        return None

# Apply the function to create a new column with the proper chain identifier.
df1['chain'] = df1['system_protein_chains_asym_id'].apply(extract_chain)

# Perform a left merge:
# - Use 'entry_pdb_id' from the first file and match it with 'PDB' from the second file.
# - Use the newly created 'chain' column from the first file and match it with 'CHAIN' from the second file.
merged_df = pd.merge(df1, 
                     df2[['PDB', 'CHAIN', 'SP_PRIMARY']], 
                     left_on=['entry_pdb_id', 'chain'], 
                     right_on=['PDB', 'CHAIN'], 
                     how='left')

# Optionally, drop the extra columns (PDB and CHAIN) from the merge.
merged_df.drop(['PDB', 'CHAIN'], axis=1, inplace=True)
# Sort the merged DataFrame by 'system_pocket_validation_average_rsrz'
merged_df.sort_values(by=['SP_PRIMARY', 'system_pocket_validation_average_rsrz'], inplace=True)


# Drop duplicates based on 'SP_PRIMARY', keeping the first occurrence (which has the lowest 'system_pocket_validation_average_rsrz' due to sorting)
merged_df = merged_df.drop_duplicates(subset='SP_PRIMARY', keep='first')

df_chembl = pd.read_csv('/Users/ak87/ownCloud/owncloud/Calculations/HACKATHON/query_results_epyc.csv')

# Save the merged output to a new CSV file.
merged_df.to_csv('/Users/ak87/ownCloud/owncloud/Calculations/HACKATHON/merged_output.csv', index=False)

# Display the resulting DataFrame.
print(merged_df)
print("hehe")

# Merge the merged_df with df_chembl on 'SP_PRIMARY' (from merged_df) and 'uniprot_id' (from df_chembl)
final_merged_df = pd.merge(merged_df, 
                           df_chembl, 
                           left_on='SP_PRIMARY', 
                           right_on='uniprot_id', 
                           how='inner')


# open coordinate csv
df_coord = pd.read_csv('/Users/ak87/ownCloud/owncloud/Calculations/HACKATHON/combined_coord.csv')

# # Merge df_coord with final_merged_df on 'system_id'
# final_merged_df = pd.merge(final_merged_df, 
#                            df_coord, 
#                            on='system_id', 
#                            how='inner')
# Save the final merged output to a new CSV file.
final_merged_df.to_csv('/Users/ak87/ownCloud/owncloud/Calculations/HACKATHON/final_merged_output.csv', index=False)

# Display the resulting DataFrame.
print(final_merged_df)

# Create a dictionary to hold the JSON structure
json_data = {}

# Group the DataFrame by SP_PRIMARY
for sp, group in final_merged_df.groupby('SP_PRIMARY'):
    print(sp)
    # Assuming that for a given SP_PRIMARY, entry_pdb_id, ligand_interacting_residues, and split are consistent
    entry_pdb_id = group['entry_pdb_id'].iloc[0]
    ligand_interacting_residues = group['ligand_interacting_residues'].iloc[0]
    # Filter df_coord for the current system_id
    system_id = group['system_id'].iloc[0]
    coord_subset = df_coord[df_coord['system_id'] == system_id]
    # print(coord_subset)
    if coord_subset.empty:
        print(f"No coordinates found for system_id: {system_id}")
        continue
    # Create a list of dictionaries for the coordinates
    coordinates = coord_subset[['resname', 'resseq', 'atom_name', 'x', 'y', 'z']].to_dict(orient='records')
    split_val = group['split'].iloc[0]
    
    # Build the list of ligands from canonical_smiles and pchembl_value
    ligands = group[['canonical_smiles', 'pchembl_value']].to_dict(orient='records')
    # Rename keys for clarity: 'canonical_smiles' becomes 'ligand_smiles' and 'pchembl_value' becomes 'affinity'
    validated_ligands = []
    for d in ligands:
        ligand_smiles = d['canonical_smiles']
        try:
            mol = Chem.MolFromSmiles(ligand_smiles)
            if mol is not None:
                # print(f"Valid SMILES: {ligand_smiles}")
                fp = AllChem.GetMorganFingerprintAsBitVect(mol,2,2048)
                validated_ligands.append({'ligand_smiles': ligand_smiles, 'fingerprint': fp, 'affinity': d['pchembl_value']})
            else:
                print(f"Invalid SMILES: {ligand_smiles}")
        except Exception as e:
            print(f"Error processing SMILES: {ligand_smiles}")
            print(e)

    # Assign the constructed dictionary to the corresponding SP_PRIMARY key
    json_data[sp] = {
        "entry_pdb_id": entry_pdb_id,
        "ligand_interacting_residues": ligand_interacting_residues,
        "binding_site_coordinates": coordinates,
        "split": split_val,
        "ligands": validated_ligands
    }

# Save the JSON object to a file
with open('/Users/ak87/ownCloud/owncloud/Calculations/HACKATHON/final_merged_output_validated_coords_fps.json', 'w') as f:
    json.dump(json_data, f, indent=4)

print("JSON file 'output.json' has been saved.")
