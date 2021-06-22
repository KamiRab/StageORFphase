#!/usr/bin/env zsh
set -euxo pipefail
arr_srr=(/media/crabier/Elements/SRA/SRR6398750.fastq /media/crabier/Elements/SRA/SRR6398751.fastq /media/crabier/Elements/SRA/SRR6398752.fastq /media/crabier/Elements/SRA/SRR6398753.fastq /media/crabier/Elements/SRA/SRR6398754.fastq /media/crabier/Elements/SRA/SRR6398755.fastq /media/crabier/Elements/SRA/SRR6398756.fastq)
for i in "${arr_srr[@]}"
do
  python ORFphase_python3.py -gff Saccer.gff -fasta Saccer3.fasta -fastq "$i" -adapt CTGTAGGCACCATCAAT
  j=$(python find_best_kmer.py -fastq "$i")
  python ORFribomap.py -gff Saccer.gff -fasta Saccer3.fasta -fastq "$i" -kmer "$j" -cutdir . -o nc -features_include nc_intergenic -options MBP
done
