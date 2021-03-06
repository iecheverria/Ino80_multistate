########################################################
# Take two fasta files, do T-coffee
# alignment and compute sequence identity.
#
# Otherwise enrich sequence and recompute alignment
# (improve coverage without reducing identity) TODO
#
# iecheverria - Sali Lab - UCSF
#########################################################

import numpy as np
import pandas as pd
import pickle
import os
import pandas as pd


def run_psiblast(fasta):
    last_it = 0
    # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
    exeline = (
        "~/SOFTW/ncbi-blast-2.7.1+/bin/psiblast -evalue 10 -gapopen 11 -gapextend 1"
    )
    exeline += " -outfmt 7 -matrix BLOSUM62 -db "
    exeline += "refseq_protein -remote -query " + fasta
    homoOut = os.popen(exeline).read()
    homoOut = homoOut.splitlines()  # now an array

    # First determine the number of iterations
    for line in homoOut:
        if "# Iteration:" in line:
            it = line.split()[2]

    # Only consider last iteration
    E = []
    for line in homoOut:
        if line == "# Iteration: " + it:
            last_it = 1
        entries = line.split()
        if len(entries) > 0:
            if entries[0] != "#" and entries[0] != "Search" and last_it == 1:
                pdb = entries[1].split("|")[3] + entries[1].split("|")[4]
                E.append(entries)
            else:
                continue
    return np.array(E)


def run_pblast(fasta1, fasta2):
    """ Run protein-protein blast"""

    names = ["evalue", "pident", "qcovs", "qstart", "qend", "sstart", "send"]
    S = pd.DataFrame(columns=names)

    exeline = "~/SOFTW/ncbi-blast-2.7.1+/bin/blastp -remote "
    exeline += "-query " + fasta2 + " -subject " + fasta1
    exeline += ' -outfmt "7 evalue pident qcovs qstart qend sstart send" '
    print(exeline)
    pair_blast = os.popen(exeline).read()
    pair_blast = pair_blast.splitlines()

    for line in pair_blast:
        print(line)
        if line[0] != "#":
            vals = line.split()
            S = S.append(pd.Series(vals, index=names), ignore_index=True)
    return S


def run_tcoffee(fasta1, fasta2, name="results"):
    os.system(f"cat {fasta1}  {fasta2} > seq_for_aln.fasta")
    os.system(
        f"t_coffee seq_for_aln.fasta -output=score_html clustalw_aln fasta_aln score_ascii phylip -maxnseq=150 -maxlen=10000 -case=upper -seqnos=off -outorder=input -run_name={name} -multi_core=4 -mode=regular > aln_{name}.txt"
    )


class clustal_to_mat:
    def __init__(self, S1_name, S2_name, file_aln, file_scores):
        self.S1_name = S1_name
        self.S2_name = S2_name
        
        self.ch1 = S1_name
        self.ch2 = S2_name

        print('NAME_______', self.S1_name, self.S2_name)
        print('File aln', file_aln)

        self.seqS1, self.seqS2 = self.read_clustal(file_aln)
        self.scores = self.read_scores(file_scores)
        print('---- seqs', self.seqS1, self.seqS2)
        
        file_out = f"mat_{S1_name}_{S2_name}.align"

        self.to_mat()
        new_mat, segments = self.identify_segments()

        file_segment = f"segments_{S1_name}_{S2_name}.csv"
        segments.to_csv(file_segment)
        self.write_files(file_out, new_mat)

    def read_clustal(self, file):
        """
        Input
        ------
        Alignment and scores files in ClustalW
        format

        Returns
        -------
        Aligment list
        """

        seqS1 = []
        seqS2 = []

        for line in open(file, "r"):
            vals = line.split()
            if len(vals) == 2 and ("*" not in vals or ":" not in vals):
                id = vals[0]
                seq = vals[1]
                if self.S1_name in id:
                    seqS1 += list(seq.split("\n")[0])
                if self.S2_name in id:
                    seqS2 += list(seq.split("\n")[0])

        return seqS1, seqS2

    def read_scores(self, file):
        """
        Get alignement scores at
        residue level

        """
        scores = []
        for line in open(file, "r"):
            vals = line.split()
            if len(vals) == 2 and vals[0] == "cons":
                scores += list(vals[1].split("\n")[0])

        return scores

    def to_mat(self):
        """
        Convert ClustalW format to 
        flat list

        """
        D1 = self.seqS1
        D2 = self.seqS2

        D = np.array([D1, D2])
 
        Sc = self.scores
        
        r1, r2 = 0, 0
        self.align = []
        self.align_high = []
        for i in range(len(D[0])):
            if D[0][i] != "-":
                r1 = r1 + 1
            if D[1][i] != "-":
                r2 = r2 + 1
                
            if D[0][i] != "-" and D[1][i] != "-":
                self.align.append(
                    [D[0][i], self.ch1, r1, D[1][i], self.ch2, r2, Sc[i]]
                )
            if (
                D[0][i] != "-"
                and D[1][i] != "-"
                and (int(Sc[i]) >= 7)
            ):
                self.align_high.append(
                    [D[0][i], self.ch1, r1, D[1][i], self.ch2, r2, Sc[i]]
                )

    def identify_segments(self, segment_length=10):
        """
        Only keep aligned segments longer
        than segment_lenght

        """
        
        names = ["qstart", "qend", "sstart", "send"]
        Segments = pd.DataFrame(columns=names)

        r1 = int(self.align_high[0][2])
        r2 = int(self.align_high[0][5])
        new_mat = []
        temp = []
        for k, entry in enumerate(self.align_high):
            if int(entry[2]) == (r1 + 1) and int(entry[5]) == (r2 + 1):
                temp.append(entry)
            else:
                if len(temp) >= segment_length:
                    new_mat.append(temp)
                    Segments = Segments.append(
                        pd.Series(
                            [temp[0][2], temp[-1][2], temp[0][5], temp[-1][6]],
                            index=names,
                        ),
                        ignore_index=True,
                    )
                    temp = []
                else:
                    temp = []
            r1 = int(entry[2])
            r2 = int(entry[5])
        # Check last segment
        if len(temp) >= 10:
            new_mat.append(temp)
            Segments = Segments.append(
                pd.Series(
                    [temp[0][2], temp[-1][2], temp[0][5], temp[-1][5]], index=names
                ),
                ignore_index=True,
            )

        # Flatten list
        new_mat = [val for sublist in new_mat for val in sublist]

        return np.array(new_mat), Segments

    def write_files(self, file_out, mat):
        np.savetxt(file_out, mat, fmt="%s")


#######################################


def main():
    prot_top = {"RUVB1":"Ruvb1", "RUVB2":"Ruvb2", "IES6":"Ies6",
                "ARP5":"Arp5", "INO80":"Ino80", "H2A":"H2A",
                "H2B":"H2B", "H3":"H3", "H4":"H4"}
    
    all_prots = ["RUVB1", "RUVB2", "IES6", "ARP5", "H2A", "H2B", "H3", "H4", "INO80"]


    prot_top = {"IES2":"IES2"}
    
    all_prots = ["IES2"]
    
    for prot in all_prots:

        fasta1 = f"/Users/iecheverria/Dropbox/UCSF/Ino80/seqs/y{prot}.fasta"
        if prot == "H2A" or prot == "H2B" or prot == "H3" or prot == "H4":
            fasta2 = f"/Users/iecheverria/Dropbox/UCSF/Ino80/seqs/h{prot}.fasta"
        else:
            fasta2 = f"/Users/iecheverria/Dropbox/UCSF/Ino80/seqs/t{prot}.fasta"

        S = run_pblast(fasta1, fasta2)
        S.to_csv(f"segments_pblast_{prot}.csv")
        run_tcoffee(fasta1, fasta2, name=prot)
        
        if prot == "H2A" or prot == "H2B" or prot == "H3" or prot == "H4":
            clustal_to_mat(
                f"{prot}_YEAST",
                f"{prot}_HUMAN",
                f"{prot}.clustalw_aln",
                f"{prot}.score_ascii",
            )
            print("--------- S blast ----------")
            print(S)

            # Covert to names in the topology file
            prot_t = prot_top[prot]
            sed_comm = f"sed -e 's/{prot}_YEAST/y{prot_t}/g' -e 's/{prot}_HUMAN/h{prot_t}/g'  mat_{prot}_YEAST_{prot}_HUMAN.align >  mat_y{prot_t}_h{prot_t}.align"
            os.system(sed_comm)

        else:
            clustal_to_mat(
                f"{prot}_YEAST",
                f"{prot}_6FML",
                f"{prot}.clustalw_aln",
                f"{prot}.score_ascii",
            )
            print("--------- S blast ----------")
            print(S)

            # Covert to names in the topology file
            prot_t = prot_top[prot]
            print(prot, prot_t)
            sed_comm = f"sed -e 's/{prot}_YEAST/y{prot_t}/g' -e 's/{prot}_6FML/t{prot_t}/g'  mat_{prot}_YEAST_{prot}_6FML.align >  mat_y{prot_t}_t{prot_t}.align"
            print(sed_comm)
            os.system(sed_comm)
        
main()
