import json
import sys

def reduce_file_size(input_file, output_file):
    with open(input_file, 'r') as f:
        data = json.load(f)
    
    # Function to convert coordinates to integers (multiply by 10 and round)
    def coords_to_ints(x, y, z):
        return [
            int(round(float(x) * 10)),
            int(round(float(y) * 10)),
            int(round(float(z) * 10))
        ]
    
    # Process all entries into a minimal format
    reduced_data = []
    if isinstance(data, list):
        for entry in data:
            if 'atoms' in entry:
                # Get only CA atoms
                ca_atoms = [atom for atom in entry['atoms'] if atom.get('atom_name') == 'CA']
                if ca_atoms:
                    # Store just the ID and integer coordinates
                    reduced_data.append({
                        'i': entry.get('system_id', '').split('__')[0],  # Keep only PDB ID
                        'c': [coords_to_ints(
                            atom.get('x', 0),
                            atom.get('y', 0),
                            atom.get('z', 0)
                        ) for atom in ca_atoms]
                    })
    else:
        if 'atoms' in data:
            ca_atoms = [atom for atom in data['atoms'] if atom.get('atom_name') == 'CA']
            if ca_atoms:
                reduced_data.append({
                    'i': data.get('system_id', '').split('__')[0],
                    'c': [coords_to_ints(
                        atom.get('x', 0),
                        atom.get('y', 0),
                        atom.get('z', 0)
                    ) for atom in ca_atoms]
                })
    
    # Write reduced file with minimal whitespace
    with open(output_file, 'w') as f:
        json.dump(reduced_data, f, separators=(',', ':'))

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python reduce_json.py input_file output_file")
        sys.exit(1)
        
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    reduce_file_size(input_file, output_file)