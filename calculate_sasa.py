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
from Bio.PDB.DSSP import DSSP


d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
  'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
  'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
  'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

# read pdb file
def get_pdb_details(pdb_files, fdir, resid, aa, atomid, cd):

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
    
    # iterate through models
    for model in structures:
      chain_list = Selection.unfold_entities(model, "C")
      structure = freesasa.Structure(pdb_files[i])
      result = freesasa.calc(structure)
      
      count = 0

      # d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
      #   'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
      #   'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
      #   'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

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
            aa_list.append(residue)

            selections = freesasa.selectArea(('%s, resn %s and resi %d and chain %s and name %s' % (resid, aa.lower(), segid, chainid, atomid),
            '%s, resn %s and resi %d and chain %s and name %s' % (resid, aa.lower(), segid, chainid, atomid)), structure, result)
      
            if resname in d.keys():
              try:
                  aa_dict[new_pdb.replace('PDB', '').upper() + '_' + chainid + '_' + d[resname] + str(segid)] = selections[resid]
              except:
                  continue

      os.chdir(fdir)
      header = 'pdb_identifier,solvent_accessibility'
      write_file(pdb.replace('.pdb', '') + '_freesasa.csv', header, aa_dict, False)
      os.chdir(cd)

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
  parser.add_argument('-fdir', '--freesasa_dir', dest='fdir', nargs='?', default="freesasa_files", type=str, help='default freesasa_files')
  parser.add_argument('-aa', '--amino_acid', dest='aa', nargs='?', default="CYS", type=str, help='default CYS, options: TYR, LYS, HIS')
  parser.add_argument('-rid', '--resid', dest='rid', nargs='?', default="cysteine", type=str, help='default cysteine, options: tyrosine, lysine, histidine')
  parser.add_argument('-aid', '--atomid', dest='aid', nargs='?', default="SG", type=str, help='default SG, options: OH, NZ, NE2')
  parser.add_argument('-wo', '--write_outfile', dest='write_outfile', nargs='?', default="True", type=str, help='default True')
  parser.add_argument('-o', '--output', dest='o', nargs='?', default="sasa_identifiers.csv", type=str, help='default sasa_identifiers.csv')
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
  output_dir = cd + '/' + args.fdir.replace('files', aa.lower() + "_files")
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

  os.chdir(output_dir)
  if args.write_outfile == 'True':
    downloaded_csv_files = [f for f in listdir(os.getcwd()) if isfile(join(os.getcwd(), f))]
    write_concat_output_file("pdbs_to_sasa.csv", downloaded_csv_files)
    write_identifier_output_file(args.o, cd)

main()
