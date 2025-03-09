# MaxBind-AI: Maximum Binding Affinity Predictor

## Overview
MaxBind-AI is a machine learning tool designed to predict the maximum achievable binding strength between small molecules and protein binding sites. This innovative approach aids drug discovery researchers in understanding the theoretical limits of binding affinity optimization for specific protein targets.

## Problem Statement
In drug development, optimizing binding affinity is crucial for:
- Reducing required drug dosages
- Minimizing potential side effects
- Improving overall drug efficacy

However, researchers often invest significant time and resources in optimizing candidate drugs without knowing the theoretical maximum binding strength achievable for a specific protein target. MaxBind-AI addresses this challenge by providing predictive insights into the maximum potential binding affinity, helping researchers make informed decisions about optimization efforts.

## Technical Approach
MaxBind-AI leverages state-of-the-art machine learning techniques to analyze protein binding site characteristics and predict maximum achievable binding affinities. The model is trained on comprehensive datasets from:
- ChEMBL: Bioactivity data from scientific literature
- PDB (Protein Data Bank): Structural information about protein-ligand complexes

## Dataset Description

### Example Data (102m__1__1.A__1.C)
- Contains example data from PLINDER database
- Features sperm whale myoglobin A with cofactor (ligand)
- Includes structural files:
  - receptor.pdb/cif: Protein structure files
  - sequences.fasta: Protein sequence data
  - chain_mapping.json: Chain identification mapping
  - water_mapping.json: Water molecule mapping
  - ligand_files/: Directory containing ligand structures

### ChEMBL Data (final_merged_output files)
- Processed ligand data from ChEMBL database
- Various preprocessing stages available
- Includes:
  - Validated coordinates
  - Training data splits
  - Raw unprocessed data

## Key Features
- Prediction of maximum achievable binding affinity for specific protein targets
- Analysis of binding site characteristics
- Integration of structural and biochemical data
- Research-guided optimization recommendations

## Data Pipeline
The project includes robust data processing pipelines for:
1. Extraction and validation of protein-ligand interaction data
2. Integration of structural information from PDB
3. Feature engineering for machine learning
4. Model training and validation

## Installation
```bash
# Clone the repository
git clone https://github.com/yourusername/MaxBind-AI.git

# Install dependencies
pip install -r requirements.txt

# Initialize the environment
./init.sh
```

## Project Structure
```
MaxBind-AI/
├── Data/                      # Processed and raw data files
│   ├── final_merged_output_validated.json        # Validated data
│   ├── final_merged_output_validated_coords.json # Data with coordinates
│   └── merged_outputs_not_processed.csv         # Raw data
├── notebooks/                 # Jupyter notebooks for analysis
│   ├── merged_explore.ipynb   # Data exploration and analysis
│   └── pdb_explore.ipynb     # PDB data analysis
├── 102m__1__1.A__1.C/        # Example protein-ligand complex
├── requirements.txt          # Python dependencies
├── reduce_json.py           # JSON processing utility
├── split_json.py           # Data splitting utility
└── init.sh                 # Environment setup script
```

