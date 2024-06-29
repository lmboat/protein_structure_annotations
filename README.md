# Protein Structure Annotations

## About
PDB structures for proteins with PDB annotations can be downloaded and parsed for amino acid numbering and residue names. SIFTS files, providing residue level mapping between PDB sequences and protein sequences, were downloaded for each PDB. Cysteines resolved in each PDB were mapped to their appropriate UniProt protein sequence and identifiers for PDB to UniProt pairs were created: PDB_C#_UniProtKBID_C#. 

## Set Up
1. Import modules
```
pip3 install biopython xmlschema freesasa
```
2. Prepare a text file with PDB identifiers <pdbs_to_download.txt>

## Usage

### Calculate the solvent accessibility of each residue in a list of PDBs according to the FreeSASA package
1. Move into the solvent accessibility directory
```
cd solvent_accessibility_calculations
```
2. Download a list of PDBs
```
python3 ../download_pdbs.py
```
3. Calculate the solvent accessibilities
```
python3 ../calculate_sasa.py
```

### Identify if each residue in a list of PDBs is involved in a disulfide bond
1. Move into the disulfide directory
```
cd disulfide_bonds
```
2. Download a list of PDBs
```
python3 ../download_pdbs.py
```
3. Identify disulfide bonds
```
python3 ../calculate_disulfides.py
```

### Map PDB residues to UniProt protein sequences
1. Move into the mapping directory
```
cd pdb_protein_mapping
```
2. Download and map SIFTS
```
python3 ../parse_sifts.py
```

## Citation
1. Boatner LM, Palafox MF, Schweppe DK, Backus KM. CysDB: a human cysteine database based on experimental quantitative chemoproteomics. Cell Chem Biol. 2023 Jun 15;30(6):683-698.e3. doi: 10.1016/j.chembiol.2023.04.004. Epub 2023 Apr 28. PMID: 37119813; PMCID: PMC10510411.
