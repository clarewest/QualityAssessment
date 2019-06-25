# coding: utf-8
import os, sys, warnings, glob, argparse 
import Bio.PDB
import numpy as np
from Bio.PDB import PDBParser, PPBuilder, PDBIO

class Target:
        def __init__(self, pdbcode, terminus, region, modeldir, datadir, reference):
            self.pdbcode = pdbcode
            self.region = region
            self.modeldir = modeldir+"/" if modeldir else self.pdbcode+"/"
            self.datadir = datadir+"/" if datadir else "./"
            self.sequence = self.get_sequence()
            self.length = len(self.sequence)
            self.seg = self.get_seg()
            ## name of file to score models against 
            ## if no reference is given, default to the name of the target
            self.pdbfile = datadir + reference+".pdb" if reference else datadir + self.pdbcode+".pdb"
            self.terminus = terminus if terminus else self.get_terminus()
            self.templen = self.seg + 15
            self.chain = "A"
            self.decoychain = "A"
            self.sampled = self.get_sampled_range()
            
        ## Read sequence from fasta file and calculate length
        def get_sequence(self):
            try:
                fastafile = self.datadir + self.pdbcode + ".fasta.txt"
                with open(fastafile, "r") as f:
                    f.readline()
                    seq = f.readline().strip()
                    return seq
            except IOError:
                print("Could not read fasta file:", fastafile)
                sys.exit()
                
        ## Read the length of the segment region (unsampled structure)
        ## unless the entire structure is being scored
        def get_seg(self):
            if self.region == "all":
                return(0)
            else:
                try:
                    segfile = self.datadir + self.pdbcode + ".seg"
                    with open(segfile, "r") as f:
                        seg = int(f.readline().strip())
                        return(seg)
                except IOError:
                    print("Could not read seg from file", segfile)
                    sys.exit()
                    
        ## Read terminus from file if not provided as an argument
        def get_terminus(self):
            termfile = self.datadir + self.pdbcode + ".term"
            try:
                with open(termfile, "r") as f:
                    return(f.readline().strip())
            except IOError:
                print("Terminus not provided as argument or via", termfile)
                sys.exit()

        ## Get the indices of residues that were sampled
        def get_sampled_range(self):
            if self.terminus == "N":
                sampled_range = range(0, self.length - self.seg)
            else:
                sampled_range = range(self.seg, self.length)
            return(sampled_range)
        
        ## Get the indices of residues to be scored
        def get_score_residues(self):
            if self.region == "all": 
                ## Score all residues
                score_residues = [True] * self.length
            elif self.region == "flex":
                ## Just score the flexible region
                if self.terminus == "N":
                    flex_range = range(self.length - self.templen, self.length - self.seg)
                else: 
                    flex_range = range(self.seg, self.templen)
                score_residues = [ True if i in flex_range else False for i in range(0, self.length)]
            else:
                ## Score only sampled residues
                score_residues = [ True if i in self.sampled else False for i in range(0, self.length)]
            return(score_residues)
        
        ## Get the indices of scaffold residues to superimpose
        def get_scaffold_residues(self):
            if self.region in [ "all", "sampled" ]:
                ## superimpose the same region as scoring
                scaffold_range = self.get_score_residues() 
            if self.terminus == "N":
                scaffold_range = range(self.length - self.seg, self.length)
            else:
                scaffold_range = range(0, self.seg)
            scaffold_residues = [ True if i in scaffold_range else False for i in range(0, self.length) ]
            return(scaffold_residues)

def set_atom_lists(ref_chain, alt_chain, seq, val):
    ref_atoms = []
    alt_atoms = []
    # Get residue numbers present in the reference structure
    resnums = [ residue.id[1] for residue in ref_chain.get_residues() if residue.id[1] <= len(seq) ]   
    # Only include equivalent residues from the model
    alt_chain_eq = [ alt_chain[i] for i in resnums ]
    # Only include equivalent from the sequence
    seq_eq = [ seq[i-1] for i in resnums ]          
    # Do the same for reference to remove eg HOH at end
    val_eq = [ val[i-1] for i in resnums ]
    ref_chain_eq = [ ref_chain[i] for i in resnums ]     
    for ref_res, alt_res, amino, allow in zip(ref_chain_eq, alt_chain_eq, seq_eq, val_eq):
        if not ref_res.resname == alt_res.resname:
            print(ref_res.resname, alt_res.resname)
            sys.exit()
        assert ref_res.id == alt_res.id
        assert amino == Bio.PDB.Polypeptide.three_to_one(ref_res.resname)
        if allow:
            for atomtype in ["CA","C","O","N"]:
                ref_atoms.append(ref_res[atomtype])
                alt_atoms.append(alt_res[atomtype])
    return ref_atoms, alt_atoms

def get_rmsd(ref_atoms, alt_atoms):
    rmsd = np.sqrt(sum([ (pair[0] - pair[1]) ** 2 for pair in zip(alt_atoms, ref_atoms) ]) / len(alt_atoms))     
    return rmsd

def score_target(target, nosuperimpose_flag):
    
    print("Scoring target:", target.pdbcode)
    print("Reference structure:", target.pdbfile)

    score_residues = target.get_score_residues()
    scaf_residues = target.get_scaffold_residues()

    print("Scoring", sum(score_residues), "residues of", target.length)
    
    io = PDBIO()
    warnings.filterwarnings('always', message='.*discontinuous at.*')
    pdb_parser = PDBParser(PERMISSIVE=1,QUIET=True)

    structure = pdb_parser.get_structure(target.pdbcode, target.pdbfile)
    nat_chain = structure[0][target.chain]

    ###########################################

    with open(target.datadir + target.pdbcode + ".lst", "r") as f:
        decoys = [ target.modeldir + "/" + line.strip() for line in f ]

    if target.region == "end":
        suffix = ".endscores"
    elif target.region == "flex":
        suffix = ".flexscores"
    elif target.region == "tempflex":
        suffix = ".tempflexscores"
    else:
        suffix = ".sampledscores"
    
    ndecoys = len(decoys)
    print("Number of decoys to score:", ndecoys)

    with open(target.pdbcode + suffix,'a') as fout:
        for count, decoy in enumerate(decoys):
            decoyname = decoy.rsplit('/')[-1]
            try:
                decoy_structure = pdb_parser.get_structure(target.pdbcode, target.datadir + decoy)
            except FileNotFoundError:
                print("ERROR File not found:", target.datadir + decoy) 
            try:
                decoy_chain = decoy_structure[0]
            except KeyError:
                print(decoy)
                continue
            decoy_chain = decoy_structure[0][target.decoychain]

            ## List of residues to score RMSD
            nat_score_atoms, decoy_score_atoms =  set_atom_lists(nat_chain, decoy_chain, target.sequence, score_residues)

            ## Superimpose if required ########
            if nosuperimpose_flag:
                segment_rmsd = 0
            else:
                ## list of atoms to set as scaffold
                nat_scaf_atoms, decoy_scaf_atoms =  set_atom_lists(nat_chain, decoy_chain, target.sequence, target.get_scaffold_residues())
                ## superimpose the scaffold of model and the reference structure
                super_imposer = Bio.PDB.Superimposer()
                super_imposer.set_atoms(nat_scaf_atoms, decoy_scaf_atoms)
                super_imposer.apply(decoy_chain.get_atoms())  
                ## get the rmsd of the scaffold against the reference segment
                ## this should be the same for all decoys
                segment_rmsd = super_imposer.rms

    ############ Get RMSD of sampled region ########
            ## Get the RMSD of the scoring region when segment is superimposed
            global_rmsd = get_rmsd(nat_score_atoms, decoy_score_atoms)

            ## Get the RMSD of the scoring region when minimised against reference
            # List of scoring atoms to set as scaffold
            nat_score_atoms, decoy_score_atoms =  set_atom_lists(nat_chain, decoy_chain, target.sequence, score_residues)
            # Superimpose scoring atoms of model and reference structure and get RMSD
            super_imposer = Bio.PDB.Superimposer()
            super_imposer.set_atoms(nat_score_atoms, decoy_score_atoms)     
            super_imposer.apply(decoy_chain.get_atoms()) 
            local_rmsd = super_imposer.rms

            ## Output RMSD of model
            fout.write("{}\t{:.3f}\t{:.3f}\t{:.3f}\n".format(decoyname, segment_rmsd, global_rmsd, local_rmsd))    

            if (ndecoys - count) % 100 == 0:
                print("Remaining:", ndecoys - count)


if __name__ == '__main__':
    ### Parse arguments ###
    ## TODO add defaults for modelpath and datapath
    parser = argparse.ArgumentParser(
            description="Score RMSDs of sampled and/or flexible regions of SAINT2 models")
    parser.add_argument("pdbcode", 
                        help="the name of the target")
    parser.add_argument("datapath", 
                        help="path to data")
    parser.add_argument("modelpath", 
                        help="the directory containing the model structures")
    parser.add_argument("--nosuperimpose",
                        help="do not perform superposition (i.e. the segment and \
                        reference structure are the same", 
                        action="store_true")
    parser.add_argument("--terminus", 
                        help="terminus of sampled region. Default: reads .term", 
                        choices = [ "C" , "N" ]) 
    parser.add_argument("--region", 
                        help="region of model to score",
                        choices = [ "flex", "sampled", "all" , "tempflex", "end"],
                        default = "flex")
    parser.add_argument("--reference", 
                        help="reference structure to score against")
    parser.add_argument("--chain", 
                        help="chain of reference structure",
                        default = "A")

    args = parser.parse_args()

    target = Target(pdbcode = args.pdbcode,
                    terminus = args.terminus,
                    region = args.region,
                    modeldir = args.modelpath + "/",
                    datadir = args.datapath + "/",
                    reference = args.reference)
    
    score_target(target, args.nosuperimpose)
