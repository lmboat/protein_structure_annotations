# Date Creadted: 230130
# Date Modified: 230216

import sys, argparse, os, subprocess, glob, re, shutil
from os.path import isfile, join, isdir
from os import listdir, path

import freesasa
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import Bio
from Bio.PDB import *
from Bio.PDB import PDBParser, NeighborSearch

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
  'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
  'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
  'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

# read pdb file
def get_pdb_details(pdb_files, ddir, resid, aa, atomid, cd):

  cd = os.getcwd()
  
  master_list = []
  count = 0

  for i in range(len(pdb_files)):
    print(pdb_files[i])
    pdb = pdb_files[i]
    new_pdb = pdb.replace('.pdb', '').upper()

    # read pdb file
    parser = Bio.PDB.PDBParser(QUIET=True) 
    try:
      structures = parser.get_structure(pdb_files[i][:-4], pdb)
    except:
      print("error in reading structure")
      continue

    sg_atoms = []
    
    # iterate through models
    for model in structures:
      chain_list = Selection.unfold_entities(model, "C")

      aa_dict = {}

      for chain in chain_list:
        chainid = chain.get_id()
        res_dict = {}
        aa_list = []
        residues = chain.get_residues()
        res_dict = {}
        chain_dict = {}

        for residue in residues:
          resname = residue.get_resname()
          segid = list(residue.id)[1]

          # skip amino acids with negative indices
          if (residue.get_resname() == aa) & (segid >= 0):

            if 'SG' in residue:
              sg_atoms.append(residue['SG'])

            # disulfide_bridge = find_close_cysteines(structures)

            if resname in d.keys():
              try:
                  aa_dict[new_pdb.replace('PDB', '').upper() + '_' + chainid + '_' + d[resname] + str(segid)] = selections[resid]
              except:
                  continue

    # Find SG atoms within 3 angstroms of each other
    ns = NeighborSearch(sg_atoms)
    disulfide_bonds = []
    for sg_atom in sg_atoms:
        neighbors = ns.search(sg_atom.coord, 3.0)
        for neighbor in neighbors:
            if neighbor != sg_atom:
                residue1 = sg_atom.get_parent()
                residue2 = neighbor.get_parent()
                chain1 = residue1.get_parent()
                chain2 = residue2.get_parent()
                bond = tuple(sorted(((chain1.get_id(), residue1.get_id()[1]), (chain2.get_id(), residue2.get_id()[1]))))
                if bond not in disulfide_bonds:
                  disulfide_bonds.append(bond)

    disulfide_ids = []
    first_cys = []
    second_cys = []

    # Output the disulfide bonds with chain information
    for bond in disulfide_bonds:
      res1 = bond[0]
      res2 = bond[1]

      if res1 != res2:
        # print(f"Disulfide bond between residue {res1[1]} in chain {res1[0]} and residue {res2[1]} in chain {res2[0]}")
        disulfide_ids.append(new_pdb + '_' + res1[0] + '_C' + str(res1[1]) + '_' + new_pdb + '_' + res2[0] + '_C' + str(res2[1]))
        first_cys.append(new_pdb + '_' + res1[0] + '_C' + str(res1[1]))
        second_cys.append(new_pdb + '_' + res2[0] + '_C' + str(res2[1]))

    disulfide_df = pd.DataFrame()
    disulfide_df['disulfide_identifier'] = disulfide_ids
    disulfide_df['first_cysteine'] = first_cys
    disulfide_df['second_cysteine'] = second_cys

    os.chdir(ddir)
    disulfide_df.to_csv(pdb.replace('.pdb', '') + '_disulfide.csv', index = False)

def list_to_string(lst, symbol):
  return (symbol.join([str(elem) for elem in lst]))

def write_file(file, header, misc, special):
  out_file = open(file, 'w')
  out_file.write(header + '\n')

  if isinstance(misc, dict):
    write_dict_file(out_file, misc)

  elif isinstance(misc, list):
    write_list_file(out_file, misc, special)

  else:
    print("Cannot write file, unrecognized format. " + str(misc))

  out_file.close()

def write_dict_file(out_file, misc):
  for k in misc:
    out_file.write(k + ',' + str(misc[k]) + '\n')

def write_list_file(out_file, misc, special):
  for i in range(len(misc)):
    current = misc[i]
    if special == True:
      st = current.output()
      out_file.write(st + '\n')
    else:
      st = ''
      for j in range(len(current)):
        st +=  str(current[j]) + ','
      out_file.write(st[:-1] + '\n')

def write_concat_output_file(output_filename, processed_sasa_files):
    # Open a new file to write the concatenated content
    with open(output_filename, 'w') as outfile:
        for fname in processed_sasa_files:
            with open(fname) as infile:
                outfile.write(infile.read())
                outfile.write("\n")  # Optionally add a newline between files

def write_identifier_output_file(output_filename, cd):
    print('Writing outfile')
    output_df = pd.read_csv("pdbs_to_sasa.csv")
    output_df = output_df.dropna()

    output_df = output_df[(output_df['pdb_identifier'] != 'pdb_identifier')]

    os.chdir(cd)
    output_df.to_csv(output_filename, index = False)

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-pdir', '--pdb_dir', dest='pdir', nargs='?', default="pdb_files", type=str, help='default pdb_files') 
  parser.add_argument('-ddir', '--disulfide_dir', dest='ddir', nargs='?', default="disulfide_files", type=str, help='default disulfide_files')
  parser.add_argument('-aa', '--amino_acid', dest='aa', nargs='?', default="CYS", type=str, help='default CYS, options: TYR, LYS, HIS')
  parser.add_argument('-rid', '--resid', dest='rid', nargs='?', default="cysteine", type=str, help='default cysteine, options: tyrosine, lysine, histidine')
  parser.add_argument('-aid', '--atomid', dest='aid', nargs='?', default="SG", type=str, help='default SG, options: OH, NZ, NE2')
  parser.add_argument('-wo', '--write_outfile', dest='write_outfile', nargs='?', default="True", type=str, help='default True')
  parser.add_argument('-o', '--output', dest='o', nargs='?', default="disulfide_identifiers.csv", type=str, help='default disulfide_identifiers.csv')
  args = parser.parse_args()

  # Move into data folder 
  os.chdir('data')

  # Set current directory
  cd = os.getcwd()

  # Check PDB directory exists
  input_dir = cd + '/' + args.pdir
  os.makedirs(input_dir, exist_ok = True)

  # Check freesasa directory exists
  if args.aa in d.keys():
    aa = d[args.aa]
  else:
    print("Not familiar with amino acid :", args.aa)
    sys.exit()

  # Check which PDB files have already been processed
  output_dir = cd + '/' + args.ddir
  os.makedirs(output_dir, exist_ok = True)
  os.chdir(output_dir)
  fs_files = [f for f in os.listdir('.') if f.endswith(".csv")]

  # Read Input CSV with column of individual PDB ids ['6MHC', '39P0']
  print('processing ' + args.aa)
  os.chdir(input_dir)
  pdb_files = [f for f in os.listdir('.') if f.endswith(".pdb")]

  # Check which PDB files need to be processed
  process_files = []
  for i in range(len(pdb_files)):
    current = pdb_files[i].replace('.pdb', '_freesasa.csv')
    if current not in fs_files:
      process_files.append(pdb_files[i])

  print("Total PDBs from input file that need to be processed: " + str(len(process_files)))

  get_pdb_details(process_files, output_dir, args.rid, args.aa, args.aid, cd)

main()
