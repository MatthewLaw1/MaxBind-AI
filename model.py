import torch
import torch.nn as nn
from transformers import RobertaModel, RobertaTokenizer
from rdkit import Chem
from rdkit.Chem import AllChem

class ObservableActivation(nn.Module):
    """A wrapper for activation functions to observe their outputs."""
    def __init__(self, activation):
        super().__init__()
        self.activation = activation
        self.outputs = None

    def forward(self, x):
        self.outputs = self.activation(x)
        return self.outputs

class SmilesEncoder(nn.Module):
    def __init__(self, embedding_dim=256):
        super().__init__()
        self.embedding_dim = embedding_dim
        
        # SMILES processing layers
        self.conv1 = nn.Conv1d(256, 512, kernel_size=3, padding=1)
        self.conv2 = nn.Conv1d(512, 512, kernel_size=3, padding=1)
        self.conv3 = nn.Conv1d(512, embedding_dim, kernel_size=3, padding=1)
        
        self.norm1 = nn.BatchNorm1d(512)
        self.norm2 = nn.BatchNorm1d(512)
        self.norm3 = nn.BatchNorm1d(embedding_dim)
        
        self.pool = nn.AdaptiveAvgPool1d(1)
        self.dropout = nn.Dropout(0.1)
        self.activation = ObservableActivation(nn.ReLU())

    def forward(self, x):
        # x shape: (batch_size, seq_len, 256)
        x = x.transpose(1, 2).to(x.device)  # Ensure tensor is on the correct device
        
        x = self.activation(self.norm1(self.conv1(x)))
        x = self.dropout(x)
        
        x = self.activation(self.norm2(self.conv2(x)))
        x = self.dropout(x)
        
        x = self.activation(self.norm3(self.conv3(x)))
        x = self.pool(x)  # (batch_size, embedding_dim, 1)
        
        return x.squeeze(-1)  # (batch_size, embedding_dim)

class BindingEmbedder(nn.Module):
    def __init__(self, feature_dim=512, embedding_dim=256):
        super().__init__()
        self.feature_dim = feature_dim
        self.embedding_dim = embedding_dim
        
        # Encoders
        self.binding_site_encoder = SmilesEncoder(embedding_dim)
        self.ligand_encoder = SmilesEncoder(embedding_dim)
        
        # Feature processing
        self.feature_proj = nn.Sequential(
            nn.Linear(feature_dim, embedding_dim),
            nn.LayerNorm(embedding_dim),
            nn.ReLU(),
            nn.Dropout(0.1)
        )
        
        # Cross-attention for binding site and ligand interaction
        self.cross_attention = nn.MultiheadAttention(
            embed_dim=embedding_dim,
            num_heads=8,
            dropout=0.1,
            batch_first=True
        )
        
        # Final fusion layers
        self.fusion_layer = nn.Sequential(
            nn.Linear(embedding_dim * 3, embedding_dim * 2),
            nn.LayerNorm(embedding_dim * 2),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(embedding_dim * 2, embedding_dim)
        )

    def forward(self, binding_site_smiles, ligand_smiles, binding_features):
        # Encode SMILES
        binding_site_emb = self.binding_site_encoder(binding_site_smiles)
        ligand_emb = self.ligand_encoder(ligand_smiles)
        
        # Process binding features
        feature_emb = self.feature_proj(binding_features)
        
        # Cross attention between binding site and ligand
        binding_site_emb = binding_site_emb.unsqueeze(1)  # Add seq_len dimension
        ligand_emb = ligand_emb.unsqueeze(1)
        attn_output, _ = self.cross_attention(
            query=ligand_emb,
            key=binding_site_emb,
            value=binding_site_emb
        )
        attn_output = attn_output.squeeze(1)
        
        # Concatenate and fuse all embeddings
        combined = torch.cat([
            binding_site_emb.squeeze(1),
            attn_output,
            feature_emb
        ], dim=-1)
        
        # Generate final embedding
        output = self.fusion_layer(combined)
        
        return output

class MaxBindModel(nn.Module):
    def __init__(self, feature_dim=512, embedding_dim=256):
        super().__init__()
        self.binding_embedder = BindingEmbedder(feature_dim, embedding_dim)
        
        # Prediction head for affinity (optional)
        self.affinity_predictor = nn.Sequential(
            nn.Linear(embedding_dim, embedding_dim // 2),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(embedding_dim // 2, 1)
        )
    
    def forward(self, binding_site_smiles, ligand_smiles, binding_features):
        # Generate informed embedding
        embedding = self.binding_embedder(
            binding_site_smiles,
            ligand_smiles,
            binding_features
        )
        
        # Predict affinity (optional)
        affinity = self.affinity_predictor(embedding)
        
        return embedding, affinity

def process_smiles(smiles):
    """Convert SMILES to molecular fingerprint tensor"""
    mol = Chem.MolFromSmiles(smiles)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    return torch.tensor(list(fp)).to(torch.device('cuda' if torch.cuda.is_available() else 'cpu'))

def prepare_batch(binding_sites, ligands, features, device):
    """Prepare a batch of data for the model"""
    # Process SMILES strings
    binding_fps = torch.stack([process_smiles(s) for s in binding_sites])
    ligand_fps = torch.stack([process_smiles(s) for s in ligands])
    
    # Convert features to tensor
    features = torch.tensor(features).to(device)
    
    return (
        binding_fps.to(device),
        ligand_fps.to(device),
        features
    )