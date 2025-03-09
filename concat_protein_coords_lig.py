import json

# Load the JSON file
with open("final_merged_output_validated_coords.json", "r") as f:
    data = json.load(f)

# Global variables to determine the maximum lengths across all parent keys
global_max_atom_name_len = 0
global_max_int_part = 0

# First pass: compute maximum atom_name length and maximum integer part length for coordinates
for parent_key, entry in data.items():
    binding_coords = entry.get("binding_site_coordinates", [])
    for atom in binding_coords:
        atom_name = atom.get("atom_name", "")
        global_max_atom_name_len = max(global_max_atom_name_len, len(atom_name))
        for coord in ["x", "y", "z"]:
            value = atom.get(coord, 0)
            int_part_length = len(str(int(abs(value))))
            global_max_int_part = max(global_max_int_part, int_part_length)

# Compute coordinate width:
# max_int_part digits + 1 (decimal point) + 3 (decimal places) + 1 (for sign)
coord_width = global_max_int_part + 1 + 3 + 1

# Temporary storage for concatenated residue strings per parent key
temp_result = {}
all_concatenated_strings = []

# Process each parent key: group atoms by residue and build concatenated strings.
for parent_key, entry in data.items():
    binding_coords = entry.get("binding_site_coordinates", [])
    residue_groups = {}
    for atom in binding_coords:
        resname = atom.get("resname", "")
        residue_groups.setdefault(resname, []).append(atom)
    
    temp_result[parent_key] = {}
    for resname, atoms in residue_groups.items():
        # Start the concatenated string with the residue name
        concatenated = resname
        for atom in atoms:
            # Pad the atom_name with trailing zeros to match the global max length
            padded_atom_name = atom.get("atom_name", "").ljust(global_max_atom_name_len, '0')
            # Format coordinates with leading zeros using the computed coord_width and 3 decimal places
            x_str = f"{atom.get('x', 0):0{coord_width}.3f}"
            y_str = f"{atom.get('y', 0):0{coord_width}.3f}"
            z_str = f"{atom.get('z', 0):0{coord_width}.3f}"
            concatenated += padded_atom_name + x_str + y_str + z_str
        temp_result[parent_key][resname] = concatenated
        all_concatenated_strings.append(concatenated)

# Determine the maximum length among all residue concatenated strings
max_concat_length = max(len(s) for s in all_concatenated_strings) if all_concatenated_strings else 0

# Build the final output including residues, ligands, and split key
final_output = {}
for parent_key, entry in data.items():
    final_output[parent_key] = {}
    
    # Process residues: pad each residue's concatenated string with trailing zeros to achieve uniform length
    final_output[parent_key]["residues"] = {}
    for resname, concat_str in temp_result[parent_key].items():
        final_output[parent_key]["residues"][resname] = concat_str.ljust(max_concat_length, '0')
    
    # Process all ligand entries: pad each ligand_smiles with trailing zeros
    final_output[parent_key]["ligands"] = []
    ligands = entry.get("ligands", [])
    for ligand in ligands:
        ligand_smiles = ligand.get("ligand_smiles", "")
        padded_smiles = ligand_smiles.ljust(max_concat_length, '0')
        affinity = ligand.get("affinity", None)
        final_output[parent_key]["ligands"].append({
            "ligand_smiles": padded_smiles,
            "affinity": affinity
        })
    
    # Include the split key from the original JSON
    final_output[parent_key]["split"] = entry.get("split", "")

# Output the final JSON result
# print(json.dumps(final_output, indent=4))
# Output the final JSON result
with open("concatenated_residue_lig_strings.json", "w") as f:
    json.dump(final_output, f, indent=4)