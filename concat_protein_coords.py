import json

# Load the JSON file
with open("final_merged_output_validated_coords.json", "r") as f:
    data = json.load(f)

# Global variables to determine the maximum lengths across all parent keys and residues
global_max_atom_name_len = 0
global_max_int_part = 0

# First pass: compute maximum atom name length and maximum integer part length for coordinates
for parent_key, entry in data.items():
    entry = json.loads(entry) if isinstance(entry, str) else entry
    binding_coords = entry.get("binding_site_coordinates", [])
    for atom in binding_coords:
        atom_name = atom.get("atom_name", "")
        # Update maximum atom_name length
        global_max_atom_name_len = max(global_max_atom_name_len, len(atom_name))
        # Check each coordinate (x, y, z)
        for coord in ["x", "y", "z"]:
            value = atom.get(coord, 0)
            int_part_length = len(str(int(abs(value))))
            global_max_int_part = max(global_max_int_part, int_part_length)

# Compute coordinate width:
# - max_int_part digits + 1 for decimal point + 3 for decimals + 1 extra slot for sign
coord_width = global_max_int_part + 1 + 3 + 1

# Intermediate storage for concatenated residue strings per parent key
temp_result = {}

# List to collect all concatenated strings to determine final uniform length
all_concatenated_strings = []

# Process each parent key separately
for parent_key, entry in data.items():
    binding_coords = entry.get("binding_site_coordinates", [])
    # Group atoms by residue name
    residue_groups = {}
    for atom in binding_coords:
        resname = atom.get("resseq", "")
        residue_groups.setdefault(resname, []).append(atom)
    
    temp_result[parent_key] = {}
    # Build concatenated strings for each residue group
    for resname, atoms in residue_groups.items():
        concatenated = atoms[0].get("resname", "") + ","
        # Start with the residue name
        # concatenated = resname
        for atom in atoms:
            # Pad atom_name with trailing zeros to the global maximum length
            padded_atom_name = atom.get("atom_name", "").ljust(global_max_atom_name_len, '0')
            # Format coordinates with leading zeros using the computed coord_width and 3 decimal places
            x_str = f"{atom.get('x', 0):0{coord_width}.3f}"
            y_str = f"{atom.get('y', 0):0{coord_width}.3f}"
            z_str = f"{atom.get('z', 0):0{coord_width}.3f}"
            concatenated += padded_atom_name + "," + x_str + "," + y_str + "," + z_str
        temp_result[parent_key][resname] = concatenated
        all_concatenated_strings.append(concatenated)

# Determine the maximum length among all residue concatenated strings
max_concat_length = max(len(s) for s in all_concatenated_strings) if all_concatenated_strings else 0

# Pad each residue's concatenated string with trailing zeros to achieve uniform length
final_output = {}
for parent_key, residues in temp_result.items():
    final_output[parent_key] = {}
    for resname, concat_str in residues.items():
        final_output[parent_key][resname] = concat_str.ljust(max_concat_length, '0')

# Output the final JSON result
with open("concatenated_residue_strings1.json", "w") as f:
    json.dump(final_output, f, indent=4)