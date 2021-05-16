#!/usr/bin/env python3
# modified 20/02/18
# From: Chris Papadopoulos
# 
# Modifications: 
# 1. There are IGORFs that are just next to the codding part but they have a 
#    STOP codon instead of overlapping to the coding. 
#    These IGORFs MUST NOT be extended and overlap with the codding
#    otherwise we generate fasta with 2 stop codons ---> ALKFPLKL*PMF*
#
""" IGR strands for InterGenic Regions.

This script defines IGRs by a left and right genomic limits
(left < right) for each chromosome. An IGR is a region
with no overlap with coding regions, both strands included.
Detection is made from GFF3 formatted file.
"""

# Library
import argparse, pickle
from collections import namedtuple
import copy
# Local
#from script_objects import Igr, Igorf
from noncoding_objects import Igr, IgorfDetector
from sequtils import reverse_complement
from constants import *

# Constants
MIN_VIEW = namedtuple("min_view",
                      ["left", "right", "left_strand", "right_strand",
                       "left_frame", "right_frame"])
# ***************************************
# ADD "rRNA_gene","tRNA_gene" pour trouver des igorfs sans rRNA et tRNA et snoRNA
CODING_FEATS = ['gene', 'pseudogene', 'ncRNA_gene', "rRNA_gene","tRNA_gene","snoRNA_gene","snRNA_gene","transposable_element_gene"]
# ***************************************
CODON_LEN = 3
# stop +
STOPS = ["TAG", "TAA", "TGA"]
# stop reverse (-)
STOPS_REV = ["TTA", "CTA", "TCA"]

FR_SH = {
    1:{
        0:1,
        1:2,
        2:0
    }
}

OPEN_KEYWD = "OPEN"
# paramètre qui peut être changer
IGORF_MINLEN = 60 #nt


def main(args):
    """ main of the program. """
    print("\n ** %sIGRs - IGORFs detection%s **\n" % (BLITUN, NA))
    # Loading fasta pickle
    print("Loading %s" % (args.fapick))
    fapick = pickle.load(open(args.fapick, "rb"))
    # Loading gff file to define coding blocks
    print("Loading %s" % (args.gff_file))
    gff = load_gff(args.gff_file, fapick, feature_type=CODING_FEATS)
    # For each seqid, detect all IRGs + associated igORFs
    igrs = {}
    print("Detecting IGRs and IGORFs in", end='\r')
    for seqid in gff:
        print("Detecting IGRs and IGORFs in %s%15s%s" % (ORIT, seqid, NA),
              end='\r')
        # Get all IGRs of the seqid
        igrs[seqid] = detect_igrs(gff[seqid], fapick["sizes"][seqid], seqid)
        # Detect orf in all frames for each IRG
        scan_frames(fapick, igrs[seqid], args.igorf_minlen)
    print("")
    write_gff(igrs)
    print_summary(igrs)
    # One particular case: when the Igorf doesnt have stop.
    # Start of chrVI for instance, check if left limit < 0 => set it to 0.



def print_summary(igrs):
    """ """
    igrs_num = 0
    keeped = 0
    removed = 0
    for seqid in igrs:
        for igr in igrs[seqid]:
            igrs_num += 1
            for igorf in igr.childs:
                if igorf.keep_me:
                    keeped += 1
                else:
                    removed += 1
    print("\nIGRs: {igrs_num}".format(**locals()))
    print("IGORFs: {tot}".format(tot=(keeped + removed)))
    print("# of keeped IGORFs: {keeped}".format(**locals()))
    print("# of removed IGORFs: {removed}".format(**locals()))

# Added by Chris Papadopoulos (21/02/2019): 
def check_the_right_stop(igorf):
    #####################################
    size_of_my_seq = len(igorf.sequence)
    o = copy.deepcopy(igorf)
    o.update_right(fapick,"OV3_same",igorf_minlen)
    size_of_my_seq_updated = len(o.sequence)
    starting_position = size_of_my_seq_updated - size_of_my_seq
    codon_before = o.sequence[starting_position-3:starting_position]
    if codon_before not in STOPS:
    #####################################
        return( igorf.update_right(fapick,"OV3_same",igorf_minlen) )
    else:
        return( igorf )

def write_gff(igrs):
    """ Write a gff file with all IGRs and IGORFs.
    IGORFs discarded by -minlen are not written in this GFF.
    """
    for seqid in igrs:
        seqid_out = "%s.gff" % (seqid)
        with open(seqid_out, "w+") as out:
            for igr in igrs[seqid]:
                gffline = igr.gffline()
                out.write("%s\n" % (gffline))
                for igorf in igr.childs:
                    if igorf.keep_me:
                        gffline = igorf.gffline()
                        out.write("%s\n" % (gffline))

        print("%s written" % (seqid_out))

def scan_frames(fapick, igrs, igorf_minlen):
    """ For each igrs of the seqid, scan the sequence to detect ORF. """

    for igr in igrs:
        seqid_size = fapick["sizes"][igr.seqid]
        seqid_seq = fapick["sequences"][igr.seqid]

        igr_seq = igr.get_sequence(fapick)
        rev_igr_seq = reverse_complement(igr_seq)

        plus_igorfs = detect_orfs(igr_seq, igr.left_genomic, "+",
                                  1, seqid_size, seqid_seq, igr)
        minus_igorfs = detect_orfs(rev_igr_seq, igr.right_genomic, "-",
                                   -1, seqid_size, seqid_seq, igr)

        all_igorfs = {
            "+":plus_igorfs,
            "-":minus_igorfs
        }
        # Define flags for location of NO_LIMITS
        nl_left = ""
        nl_right = ""

        for strand in all_igorfs:
            for frame in all_igorfs[strand]:
                # igorf of the current strand and frame
                sf_igorfs = all_igorfs[strand][frame]
                for idx, igorf in enumerate(sf_igorfs):
                    # Flag the size parameter
                    if igorf.length < igorf_minlen:
                        igorf.keep_me = False
                    # First IGORF
                    if idx == 0:
                        if strand == "+" and igr.left_frame >= 0:
                            # If igr.left_frame<0 ==> No igr.left (startofchrom)
                            # strand +, left
                            if igr.left_strand == strand:
                                if igr.left_frame == igorf.absolute_frame:
                                    igorf.location = "UTR3"
                                else:
                                    # DO NOT CHANGE!!!
                                    igorf.update_left(fapick,
                                                      "OV3_same",
                                                      igorf_minlen)
                            else:
                                # DO NOT CHANGE
                                igorf.update_left(fapick,
                                                  "OV5_other",
                                                  igorf_minlen)
                                    
                        elif strand == "-" and igr.right_frame >= 0:
                            # If igr.right_frame<0 ==> No igr.right (endofchrom)
                            # strand -, right
                            if igr.right_strand == strand:
                                if igr.right_frame == igorf.absolute_frame:
                                    igorf.location = "UTR3"
                                else:
                                    check_the_right_stop(igorf)

                                    
                            else:
                                igorf.update_right(fapick,
                                                   "OV5_other",
                                                   igorf_minlen)
                        if len(sf_igorfs) == 1:
                            if strand == "+":
                                nl_left = igorf.location
                            else:
                                nl_right = igorf.location
                    # Last IGORF
                    if idx == len(sf_igorfs) - 1:
                        if strand == "+" and igr.right_frame >= 0:
                            # strand +, right
                            if igr.right_strand == strand:
                                if igr.right_frame == igorf.absolute_frame:
                                    igorf.location = "UTR5"
                                else:
                                    # There are cases which are next to the 
                                    # coding but have a STOP codon at the 
                                    # end of their sequence. So we will not 
                                    # search for overlapping.
                                    if igorf.sequence[-3:] not in STOPS:
                                        igorf.update_right(fapick,
                                                       "OV5_same",
                                                       igorf_minlen)
                            else:
                                # There are cases which are next to the 
                                # coding but have a STOP codon at the 
                                # end of their sequence. So we will not 
                                # search for overlapping.
                                if igorf.sequence[-3:] not in STOPS:
                                    igorf.update_right(fapick,
                                                   "OV3_other",
                                                   igorf_minlen)
                        elif strand == "-" and igr.left_frame >= 0:
                            # strand -, left
                            if igr.left_strand == strand:
                                if igr.left_frame == igorf.absolute_frame:
                                    igorf.location = "UTR5"
                                else:
                                    # There are cases which are next to the 
                                    # coding but have a STOP codon at the 
                                    # end of their sequence. So we will not 
                                    # search for overlapping.
                                    if igorf.sequence[-3:] not in STOPS:
                                        igorf.update_left(fapick,
                                                      "OV5_same",
                                                      igorf_minlen)
                                    
                            else:
                                if igorf.sequence[-3:] not in STOPS:
                                    igorf.update_left(fapick,
                                                  "OV3_other",
                                                  igorf_minlen)
                        # If len == 1, the left or right limit is already
                        # check wich one is missing
                        # and create the hybride location
                        if len(sf_igorfs) == 1:
                            if nl_left:
                                location = "%s::%s-%s" % (OPEN_KEYWD,
                                                          nl_left,
                                                          igorf.location)
                            else:
                                location = "%s::%s-%s" % (OPEN_KEYWD,
                                                          igorf.location,
                                                          nl_right)
                            nl_left = ""
                            nl_right = ""
                            igorf.location = location
                    igr.childs.append(igorf)

def detect_orfs(sequence, first_pos, strand,
                increment, seqid_size, seqid_seq, igr):
    """ Detect ORFs in the 3 frames of the sequence.
    ORF is defined as space between two stops codons.
    So, in a frame, the ORF starts at the 1st base until the first stop.
    strand is -1 ('-') or 1 ('+')
    """
    igorfs = {
        0:[],
        1:[],
        2:[]
    }

    # Exclude too small IGORFs
    if igr.length < IGORF_MINLEN:
        return igorfs

    for frame in range(3):
        frame_seq = sequence[frame:]
        frame_codons = [frame_seq[i:i+3] for i in range(0, len(frame_seq), 3)]
        # To stop at the last codon in the ORF frame
        if len(frame_codons[-1]) == CODON_LEN:
            final_idx = -1
        else:
            final_idx = -2

        tmp_orf = ""
        tmp_pos = []
        for idx, codon in enumerate(frame_codons):
            if codon not in STOPS:
                tmp_orf += codon
                tmp_pos.append(idx)
                if idx == len(frame_codons) + final_idx:
                    orf_start = first_pos\
                                + (frame * increment)\
                                + ((tmp_pos[0] * 3) * increment)
                    orf_stop = orf_start \
                               + (len(tmp_orf) * increment) \
                               - (1 * increment)
                    if increment == 1:
                        orf_absolute_frame = orf_start % 3
                    else:
                        orf_absolute_frame = (seqid_size - (orf_start - 1)) % 3
                    igorf = IgorfDetector(igr, orf_start, orf_stop, strand,
                                  increment, tmp_orf, orf_absolute_frame)
                    igorfs[orf_absolute_frame].append(igorf)
                    # End of the loop
            else:
                tmp_orf += codon
                tmp_pos.append(idx)
                orf_start = first_pos \
                            + (frame * increment) \
                            + ((tmp_pos[0] * 3) * increment)
                orf_stop = orf_start \
                           + (len(tmp_orf) * increment) \
                           - (1 * increment)
                if increment == 1:
                    orf_absolute_frame = orf_start % 3
                else:
                    orf_absolute_frame = (seqid_size - (orf_start - 1)) % 3
                igorf = IgorfDetector(igr, orf_start, orf_stop, strand,
                              increment, tmp_orf, orf_absolute_frame)
                igorfs[orf_absolute_frame].append(igorf)
                # reset
                tmp_orf = ""
                tmp_pos = []
    return igorfs

def detect_igrs(tuple_views, seqid_size, seqid):
    """ From tuple_views of coding regions in a seqid,
    this function detect all intervals between them.
    """
    # Check overlaps between coding regions without taking strand in account
    coding_blocks = check_overlaps(sorted(tuple_views))
    #print(sorted(coding_blocks))
    # IGR are the blank spaces between coding blocks.
    igrs = []

    # Loop on coding_blocks to get intern IRGs + 2 special
    # cases for extern IGRs (5', 3' of chromosome)
    for idx, cob in enumerate(coding_blocks):
        if idx < len(coding_blocks) - 1:

            if idx == 0: # right most IGR
                igr = Igr(
                    seqid,
                    1,
                    cob.left - 1,
                    left_strand=-10, #No left strand, 5' of the chromosome
                    right_strand=cob.left_strand,
                    left_frame=-10,
                    right_frame=cob.left_frame
                )
                igrs.append(igr)

            # Internal IGRs
            next_cob = coding_blocks[idx + 1]
            igr_start = cob.right + 1
            igr_end = next_cob.left - 1
            igr_lim = (igr_start, igr_end)
            igr = Igr(
                seqid,
                igr_start,
                igr_end,
                left_strand=cob.right_strand,
                right_strand=next_cob.left_strand,
                left_frame=cob.right_frame,
                right_frame=next_cob.left_frame
            )
        else:
            # Last coding region so left most igr
            igr = Igr(
                seqid,
                cob.right + 1,
                seqid_size,
                left_strand=cob.right_strand,
                right_strand=-10,  #No right strand, 3' of the chromosome,
                left_frame=cob.right_frame,
                right_frame=-10
            )
        igrs.append(igr)
    return igrs
# why recursif
def check_overlaps(tuples):
    """ Return the wider tuple if two tuples are overlapping eachother.
    Also works for included tuples.
    ex1: [(10, 20), (15, 30)] => [(10, 30)]
    ex2: [(10, 40), (20, 30)] => [(10, 40)]
    """
    for i, cob1 in enumerate(tuples):
        for j, cob2 in enumerate(tuples[i+1:], i+1):
            if cob1.right >= cob2.left:
                if cob2.right >= cob1.right: #Overlapped case
                    tuples[i] = MIN_VIEW(cob1.left, tuples.pop(j).right,
                                         cob1.left_strand, cob2.right_strand,
                                         cob1.left_frame, cob2.right_frame)
                else: #Included case
                    tuples[i] = MIN_VIEW(cob1.left, tuples.pop(i).right,
                                         cob1.left_strand, cob1.right_strand,
                                         cob1.left_frame, cob1.right_frame)
                return check_overlaps(tuples)
    return tuples


def load_gff(gff_file, fapick, feature_type=["gene"]):
    """ Load all gff file features with a minimal representation as tuple:
    (left_limit, right_limit, left_strand,
    right_strand, left_frame, right_frame)
    """
    total_features = 0
    gff = {}
    print(CODING_FEATS)
    with open(gff_file, 'r') as gff_input:
        for line in gff_input:
            if not line.startswith(GFF_COM) and not line.isspace():
                splitted_line = line.strip().split(GFF_SEP)
                seqid = splitted_line[0]
                ftype = splitted_line[2]
                left = int(splitted_line[3])
                right = int(splitted_line[4])
                strand = splitted_line[6]
                seqid_size = fapick["sizes"][seqid]


                if strand == "+":
                    # Already in 1 based index,
                    # adjust the stop to the first nt of stop codon
                    # frame recover
                    frame_left = left % 3
                    frame_right = (right - 2) % 3
                # strand - 
                else:
                    # -1 to adjust to the 1 based index,
                    # and then +1 for the stop (first nt of the stop codon)
                    frame_left = (seqid_size - (left + 1)) % 3
                    frame_right = (seqid_size - (right - 1)) % 3
                if ftype in CODING_FEATS:
                    # Strand is two times, as future coding blocks can have
                    # different strand for both extremities
                    tuple_view = MIN_VIEW(left, right,
                                          strand, strand,
                                          frame_left, frame_right)
                    total_features += 1
                    try:
                        gff[seqid].append(tuple_view)
                    except KeyError: # seqid not init in gff
                        gff[seqid] = [tuple_view]
    # gff = dictionary with list of coding regions for each seqid
    print("    > # of chromosomes loaded: %s"  % (len(gff)))
    print("    > # of features loaded: %s" % (total_features))
    return gff

def return_first(iterable):
        """ Return first value in an iterable. """
        return iterable[0]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "gff_file",
        help="Gff3 formatted file"
    )
    parser.add_argument(
        "-fapick", "--fapick",
        required=True,
        help="Python3 pickle with fasta file information."
        "Generated from fasta2pickle.py script."
    )
    parser.add_argument(
        "-minlen", "--igorf_minlen",
        type=int,
        default=60,
        help="IGORFs with length < -minlen will be discarded."
        "Default = %s" % (IGORF_MINLEN)
    )
    args = parser.parse_args()
    main(args)
# Run program
# mkdir igorf_norRNA_tRNA_snoRNA_snRNA_gene
# cd igorf_norRNA_tRNA_snoRNA_snRNA_gene
# python3.5 /mnt/workspace/Maxime_Stage_M2/software/strali/detect_igr.py /mnt/workspace/Maxime_Stage_M2/data/data_ribosome_profiling_resultat_JN/gene_de_novo/Saccer.gff -fapick /mnt/workspace/Maxime_Stage_M2/data/data_ribosome_profiling_resultat_JN/gene_de_novo/SACCER3.fapick -minlen 60 
