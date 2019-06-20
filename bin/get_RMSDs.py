# coding: utf-8
import glob
import Bio.PDB
import numpy
from Bio.PDB import PDBParser, PPBuilder, PDBIO
import warnings
import sys
import os

def set_atom_lists(ref_chain, alt_chain, seq, val):
  ref_atoms = []
  alt_atoms = []
  resnums = [ residue.id[1] for residue in ref_chain.get_residues() if residue.id[1] <= len(seq) ]   # Get residue numbers present in native
  alt_chain_eq = [ alt_chain[i] for i in resnums ]                                                   # Only include equivalent from decoy
  seq_eq = [ seq[i-1] for i in resnums ]                                                             # Only include equivalent from sequence
  val_eq = [ val[i-1] for i in resnums ]
  ref_chain_eq = [ ref_chain[i] for i in resnums ]                                                   # Do the same for native to remove eg HOH at end
  for ref_res, alt_res, amino, allow in zip(ref_chain_eq, alt_chain_eq, seq_eq, val_eq):
      if not ref_res.resname == alt_res.resname:
        print ref_res.resname, alt_res.resname
        exit()
      assert ref_res.id      == alt_res.id
      assert amino == Bio.PDB.Polypeptide.three_to_one(ref_res.resname)
      if allow:
        for atomtype in ["CA","C","O","N"]:
          ref_atoms.append(ref_res[atomtype])
          alt_atoms.append(alt_res[atomtype])
  return ref_atoms, alt_atoms

def get_rmsd(ref_atoms, alt_atoms):
    rmsd = numpy.sqrt(sum([ (pair[0] - pair[1]) ** 2 for pair in zip(alt_atoms, ref_atoms) ]) / len(alt_atoms))     
    return rmsd

myhost = os.uname()[1]
args = sys.argv
homo = "--homo" in args or "-h" in args         # Segment and native not from same file 
crystal = "--crystal" in args or "-c" in args   # Segment from a different crystal structure
end = "--end" in args or "-e" in args           # Ignore segment in RMSD
saulo = "--saulo" in args or "-s" in args       # No segment, different file set-up
rosetta = "--rosetta" in args or "-r" in args   # Saulo but with different file set-up
flex = "--flex" in args or "-f" in args         # Score the flexible region
tempflex = "--tempflex" in args or "-t" in args # Score the flexible region against the template
if crystal:
  args.remove("-c")
if homo:
  args.remove("-h")
if saulo:
  args.remove("-s")
if end:
  args.remove("-e")
if rosetta:
  args.remove("-r")
if flex:
  args.remove("-f")
if tempflex:
  args.remove("-t")
if len(args) != 4:
  print "USAGE:    script.py [-h] [-c] [-e] [-s] <pdbcode> <pathtodata> <decoysuperfolder>"
pdbcode = args[1]                                                                       # Full name of target 
pathtodata = args[2]                                                                    # Proteinprep directory
decoysuperfolder = args[3]                                                              # Decoy directory path 
rev = "_rev" == decoysuperfolder[-4:]

print "Homology model scaffold:", homo

io = PDBIO()
warnings.filterwarnings('always', message='.*discontinuous at.*')
pdb_parser = PDBParser(PERMISSIVE=1,QUIET=True)

if saulo or rosetta:
  seg = 0                                                       
else:
  with open(pathtodata + "/" + pdbcode + ".seg") as fin:
    seg = int(fin.readline().strip())
if flex or tempflex:
  with open(pathtodata + "/"+ pdbcode + ".templen") as fin:
    templen = int(fin.readline().strip())
if tempflex:
  pdb = pathtodata + "/" + pdbcode + ".pdb"
elif homo:
  pdb = pathtodata + "/real_" + pdbcode + ".pdb"                                            # PDB file of native structure (pdbcode.pdb is the template for seg)
else:
  pdb = pathtodata + "/" +  pdbcode + ".pdb"
native = pdbcode                                      
print "Using target structure: ", pdb
if tempflex:
#  with open("../../" + pdbcode + ".tempchain") as fin:
#    chain = fin.readline().strip()
  with open(pathtodata + "/" + pdbcode + ".chain") as fin:    # Template chain ID
      chain = fin.readline().strip()
      if len(chain) is 0:
        chain = " "
else:
  chain = "A"
if saulo or rosetta:                                           # Scaffold chain is A for homology model or normal decoy
  decoychain = "A"
elif crystal or tempflex:
  with open(pathtodata + "/" + pdbcode + ".chain") as fin:    # Otherwise its template crystal structure chain
      decoychain = fin.readline().strip()
      if len(decoychain) is 0:
        decoychain = " "
else:
  decoychain = chain                                        # For normal SAINT2 Eleanor it's the same as seg
#  tmhs = [ [ int(i) for i in line.strip().split() ] for line in lines[4:] ]
with open(pathtodata + "/" + pdbcode + ".length") as fin:
  length = int(fin.readline().strip())

print "Number of residues to score:", length-seg                                      # Length of sampled region 

with open(pathtodata + "/" + pdbcode + ".fasta.txt") as fin:
  lines = fin.readlines()
  sequence = lines[1].strip()

print(len(sequence))
print(length)
print(seg)

########## Setting atoms to score ##########
score_residues = [False] * len(sequence)
if saulo or rosetta:
  score_residues = [True] * len(sequence)           # Inclue all in scoring
else:
  if rev:
      samplerange = range(0, length-seg)
  else:
      samplerange = range(seg,length)
  for i in samplerange:
    score_residues[i] = True                         # Include only sampled region in scoring

######### Settings atoms as scaffold ######
if homo:                                            # If seg is not from the native, it needs to be superimposed
  homo_scaf_residues = [False] * len(sequence)
  if rev:
      scafrange = range(length-seg, length)
  else:
      scafrange = range(0,seg)
  for i in scafrange:
    homo_scaf_residues[i] = True
if end or saulo or rosetta:
  homo_scaf_residues =  [ i for i in score_residues ]  # Superimpose same region as scoring 

######## Set flex region for scoring ######            # Score just the flex region
if flex or tempflex:
  flex_score_residues = [False] * len(sequence)
  if rev:
      flexrange = range(length-templen,length-seg)
  else:
      flexrange = range(seg,templen)
  for i in flexrange:
      flex_score_residues[i] = True

###########################################

structure = pdb_parser.get_structure(native, pdb)
nat_chain = structure[0][chain]

if saulo:
  decoys = [ decoyfile for decoyfile in glob.glob(decoysuperfolder + "/" + pdbcode + "*linear*") ]
elif rosetta: 
  decoys = [ decoyfile for decoyfile in glob.glob(decoysuperfolder + "/*") ]
else:
    decoys = [ decoyfile for decoyfile in glob.glob(decoysuperfolder + "/" + myhost + "/" + pdbcode + "*linear*") if decoyfile[-4:] != ".cut" and decoyfile[-4:] != ".end" and decoyfile[-9:] != ".rescored" ]

if end:
  pref = "endscore_"
elif flex:
  pref = "flexscore_"
elif tempflex:
  pref = "tempflex_"
else:
  pref = "pyscore_"

#if os.path.exists(pref + decoysuperfolder):
#  with open(pref + decoysuperfolder) as fin:
#    donedecoys = [ line.strip().split()[0] for line in fin.readlines() ]
#  donedecoys = [ decoysuperfolder + "/" + filename.rsplit("_",2)[1] + "/" + filename for filename in donedecoys ]
#  decoys = sorted(list(set(decoys) - set(donedecoys)))


count = len(decoys)
print "Number of decoys to score:", count  
with open(pref + decoysuperfolder,'a') as fout:
  for decoy in decoys:
    decoyname = decoy.rsplit('/',1)[1]
    decoy_structure = pdb_parser.get_structure(pdbcode, decoy)
    try:
      decoy_chain = decoy_structure[0]
    except KeyError:
      print decoy
      continue
    decoy_chain = decoy_structure[0][decoychain]

    nat_score_atoms, decoy_score_atoms =  set_atom_lists(nat_chain, decoy_chain, sequence, score_residues)     # List of residues to score RMSD

############### Superimpose is required ########
    if homo or saulo or end or rosetta or tempflex:

      nat_scaf_atoms, decoy_scaf_atoms =  set_atom_lists(nat_chain, decoy_chain, sequence, homo_scaf_residues)     # List of atoms to set as scaffold

      super_imposer = Bio.PDB.Superimposer()
      super_imposer.set_atoms(nat_scaf_atoms, decoy_scaf_atoms)                                                    # Superimpose scaffold of decoy and target
      super_imposer.apply(decoy_chain.get_atoms())  
      if count == 1 and (not tempflex) and (not end):
        with open("homo_seg_accuracy",'a') as foutseg:
            foutseg.write(decoysuperfolder + " " + pref[:-1] + "\t%0.2f" % (super_imposer.rms) + "\n")                      # Output RMSD of homology model against target segment (only for the last decoy)

      #io=Bio.PDB.PDBIO()
      #io.set_structure(decoy_structure)
      #io.save("pdb_out_filename.pdb")

############ Get RMSD of sampled region ########
    if not tempflex:
      rmsd = get_rmsd(nat_score_atoms, decoy_score_atoms)                                                   # Get RMSD of sampled region after superimposition
    else:
      rmsd = super_imposer.rms
    if not (flex or tempflex):
      fout.write("{}\t{:.6f}\n".format(decoyname, rmsd))                                                                                     # Output RMSD of decoy
    else:
      nat_flex_score_atoms, decoy_flex_score_atoms =  set_atom_lists(nat_chain, decoy_chain, sequence, flex_score_residues)     # List of atoms to set as scaffold
      flex1_rmsd = get_rmsd(nat_flex_score_atoms, decoy_flex_score_atoms)
      super_imposer = Bio.PDB.Superimposer()
      super_imposer.set_atoms(nat_flex_score_atoms, decoy_flex_score_atoms)                                                    # Superimpose scaffold of decoy and target
      super_imposer.apply(decoy_chain.get_atoms()) 
      flex2_rmsd = super_imposer.rms
      if flex:
        fout.write("{}\t{:.6f}\t{:.6f}\t{:.6f}\n".format(decoyname, rmsd, flex1_rmsd, flex2_rmsd))            # Output RMSD of decoy
      elif tempflex:
        fout.write("{}\t{:.6f}\t{:.6f}\t{:.6f}\n".format(decoyname, rmsd, flex1_rmsd, flex2_rmsd))            # Output RMSD of decoy

    count -= 1
    if count % 200 == 0:
      print "count", count
  #if homo:
  #  print '\t'.join([str(res.id[1]) for res,val in zip(nat_chain,homo_valid) if val])
  #print '\t'.join([str(res.id[1]) for res,val in zip(nat_chain,valid) if val])
