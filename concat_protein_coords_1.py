import json

def convert_string_to_numbers(s):
    """
    Convert a string into a list of numbers based on each letter's position in the alphabet.
    For example, 'TYR' -> [20, 25, 18] (T=20, Y=25, R=18).
    Non-alphabet characters are ignored.
    """
    return [int(''.join([str(ord(char.upper()) - 64) for char in s if char.isalpha()]))]

def process_data(data):
    """
    Process the JSON data to group binding_site_coordinates by residue sequence (resseq).
    For each residue, create a list that starts with the numeric conversion of the resname.
    Then, for every atom in that residue, append the numeric conversion of the atom_name
    followed by the coordinates (x, y, z).
    """
    result = {}
    # Process each top-level key (e.g. 'A0A0H2UPP7')
    for parent_key, content in data.items():
        residues = {}
        for coord in content.get('binding_site_coordinates', []):
            # Use resseq (converted to string) as the grouping key
            resseq = str(coord['resseq'])
            # If this residue hasn't been processed yet, initialize its list with resname conversion
            if resseq not in residues:
                residues[resseq] = convert_string_to_numbers(coord['resname'])
            # Convert atom_name to numbers
            atom_nums = convert_string_to_numbers(coord['atom_name'])
            # Append the atom's numeric representation and its coordinates
            residues[resseq].extend(atom_nums)
            residues[resseq].extend([coord['x'], coord['y'], coord['z']])
        result[parent_key] = residues
    return result

def pad_lists(result):
    """
    Determine the maximum list length across all residues for all top-level keys.
    Pad each list with trailing zeros so that all lists have the same length.
    """
    max_len = 0
    # First, compute the maximum length
    for parent_key in result:
        for resseq, lst in result[parent_key].items():
            if len(lst) > max_len:
                max_len = len(lst)
    # Now, pad each list with zeros until its length equals max_len
    for parent_key in result:
        for resseq, lst in result[parent_key].items():
            lst.extend([0] * (max_len - len(lst)))
    return result

if __name__ == '__main__':
    # Load the input JSON file (assumed to be named 'input.json')
    with open('final_merged_output_validated_coords.json', 'r') as infile:
        data = json.load(infile)
    
    # Process the data: group by residue and convert string fields to numbers
    processed_data = process_data(data)
    
    # Pad the resulting lists with trailing zeros so that all lists have equal length
    padded_data = pad_lists(processed_data)
    
    # Write the output JSON to a file (named 'output.json') with an organized structure
    with open('output.json', 'w') as outfile:
        json.dump(padded_data, outfile, indent=4)