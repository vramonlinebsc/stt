"""
BMRB Relaxation Data Fetcher
Downloads and parses NMR relaxation data from BMRB database
"""

import requests
import pandas as pd
from pathlib import Path
import json

class BMRBFetcher:
    """Fetch and parse BMRB relaxation data"""
    
    BASE_URL = "https://bmrb.io/ftp/pub/bmrb/entry_directories/"
    API_URL = "https://api.bmrb.io/v2"
    
    def __init__(self, cache_dir="./bmrb_data"):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        
    def search_relaxation_entries(self, min_r1=True, min_r2=True, min_noe=True):
        """
        Search BMRB for entries with relaxation data
        
        Returns list of BMRB IDs with requested relaxation types
        """
        # Known good entries from literature
        well_studied_proteins = {
            'GB3': [15477, 17769],  # GB3 - extensively studied
            'ubiquitin': [6457, 15410, 19684],  # Ubiquitin
            'TDP43': [26823],  # TDP-43 C-terminal domain
        }
        
        print("Starting with well-characterized proteins:")
        for protein, ids in well_studied_proteins.items():
            print(f"  {protein}: BMRB IDs {ids}")
        
        return well_studied_proteins
    
    def fetch_entry(self, bmrb_id):
        """
        Fetch complete BMRB entry
        """
        cache_file = self.cache_dir / f"bmrb_{bmrb_id}.json"
        
        if cache_file.exists():
            print(f"Loading cached entry {bmrb_id}")
            with open(cache_file) as f:
                return json.load(f)
        
        print(f"Fetching BMRB entry {bmrb_id}...")
        url = f"{self.API_URL}/entry/{bmrb_id}"
        response = requests.get(url)
        
        if response.status_code == 200:
            data = response.json()
            with open(cache_file, 'w') as f:
                json.dump(data, f, indent=2)
            return data
        else:
            print(f"Failed to fetch {bmrb_id}: {response.status_code}")
            return None
    
    def parse_relaxation_data(self, entry_data, bmrb_id):
        """
        Extract relaxation rates from BMRB entry
        """
        if not entry_data:
            return None
        
        relaxation_data = {
            'bmrb_id': bmrb_id,
            'R1': [],
            'R2': [],
            'NOE': [],
            'CCR': []
        }
        
        # Parse different relaxation types
        # BMRB structure varies, so we need to handle multiple formats
        try:
            # Look for T1 (R1) data
            if 'heteronucl_T1_relaxation' in entry_data:
                for t1_loop in entry_data['heteronucl_T1_relaxation']:
                    if 'data' in t1_loop:
                        for row in t1_loop['data']:
                            relaxation_data['R1'].append({
                                'residue': row.get('Comp_index_ID'),
                                'atom': row.get('Atom_ID', 'N'),
                                'value': float(row.get('Val', 0)),
                                'error': float(row.get('Val_err', 0))
                            })
            
            # Look for T2 (R2) data
            if 'heteronucl_T2_relaxation' in entry_data:
                for t2_loop in entry_data['heteronucl_T2_relaxation']:
                    if 'data' in t2_loop:
                        for row in t2_loop['data']:
                            relaxation_data['R2'].append({
                                'residue': row.get('Comp_index_ID'),
                                'atom': row.get('Atom_ID', 'N'),
                                'value': float(row.get('Val', 0)),
                                'error': float(row.get('Val_err', 0))
                            })
            
            # Look for NOE data
            if 'heteronucl_NOEs' in entry_data:
                for noe_loop in entry_data['heteronucl_NOEs']:
                    if 'data' in noe_loop:
                        for row in noe_loop['data']:
                            relaxation_data['NOE'].append({
                                'residue': row.get('Comp_index_ID'),
                                'atom': row.get('Atom_ID', 'N'),
                                'value': float(row.get('Val', 0)),
                                'error': float(row.get('Val_err', 0))
                            })
            
            # Look for cross-correlation data (rarer)
            if 'cross_correlation_DD_CSA' in entry_data:
                for ccr_loop in entry_data['cross_correlation_DD_CSA']:
                    if 'data' in ccr_loop:
                        for row in ccr_loop['data']:
                            relaxation_data['CCR'].append({
                                'residue': row.get('Comp_index_ID'),
                                'type': row.get('Rex_type', 'DD_CSA'),
                                'value': float(row.get('Val', 0)),
                                'error': float(row.get('Val_err', 0))
                            })
        
        except Exception as e:
            print(f"Error parsing entry {bmrb_id}: {e}")
            return None
        
        # Summary
        print(f"\nEntry {bmrb_id} relaxation data:")
        print(f"  R1: {len(relaxation_data['R1'])} measurements")
        print(f"  R2: {len(relaxation_data['R2'])} measurements")
        print(f"  NOE: {len(relaxation_data['NOE'])} measurements")
        print(f"  CCR: {len(relaxation_data['CCR'])} measurements")
        
        return relaxation_data
    
    def get_pdb_structure(self, entry_data):
        """Extract PDB ID if available"""
        try:
            if 'related_entries' in entry_data:
                for entry in entry_data['related_entries']:
                    if entry.get('Database_name') == 'PDB':
                        return entry.get('Database_accession_code')
        except:
            pass
        return None
    
    def download_pdb(self, pdb_id, output_dir=None):
        """Download PDB structure file"""
        if output_dir is None:
            output_dir = self.cache_dir
        else:
            output_dir = Path(output_dir)
            output_dir.mkdir(exist_ok=True)
        
        pdb_file = output_dir / f"{pdb_id}.pdb"
        if pdb_file.exists():
            print(f"PDB {pdb_id} already cached")
            return str(pdb_file)
        
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(url)
        
        if response.status_code == 200:
            with open(pdb_file, 'w') as f:
                f.write(response.text)
            print(f"Downloaded PDB {pdb_id}")
            return str(pdb_file)
        else:
            print(f"Failed to download PDB {pdb_id}")
            return None

# Example usage
if __name__ == "__main__":
    fetcher = BMRBFetcher()
    
    # Get known good entries
    proteins = fetcher.search_relaxation_entries()
    
    # Fetch GB3 as test case (most studied protein in NMR)
    print("\n" + "="*60)
    print("Fetching GB3 data (BMRB 15477)...")
    print("="*60)
    
    entry = fetcher.fetch_entry(15477)
    
    if entry:
        # Parse relaxation data
        relax_data = fetcher.parse_relaxation_data(entry, 15477)
        
        # Get associated PDB structure
        pdb_id = fetcher.get_pdb_structure(entry)
        if pdb_id:
            print(f"\nAssociated PDB structure: {pdb_id}")
            pdb_file = fetcher.download_pdb(pdb_id)
            print(f"Structure saved to: {pdb_file}")
    
    print("\n" + "="*60)
    print("Setup complete! Data cached in ./bmrb_data/")
    print("="*60)
    print("\nNext steps:")
    print("1. Validate data parsing")
    print("2. Build JAX forward model")
    print("3. Implement differentiable relaxation engine")
