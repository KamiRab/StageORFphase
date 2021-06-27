#!/usr/bin/env bash
set -euxo pipefail
#Get args
kmer1=26
kmer2=30
for i in "$@"; do
  case $i in
    -gff)
      shift
      if [[ -f "$1" ]]; then
        GFF="$1"
      else
        echo "File does not exist"
      fi
      shift
      ;;
    -fasta)
      shift
      if [[ -f "$1" ]]; then
        FASTA="$1"
      else
        echo "File does not exist"
        exit 1
      fi
      shift
      ;;
    -fastq)
      shift
      if [[ -f "$1" ]]; then
        FASTQ="$1"
      else
        echo "File does not exist"
        exit 1
      fi
      shift
      ;;
    -adapt)
      shift
        ADAPT="$1"
      shift
      ;;
    -kmer)
      shift
      if [[ "$1" =~ ^[0-9]+$ ]];then
        kmer1="$1"
        exit 1
      else
        echo "Not a number"
      fi
      shift
      if [[ "$1" =~ ^[0-9]+$ ]];then
        kmer2="$1"
      else
        echo "Not a number"
        exit 1
      fi
      shift
      ;;
    -h|--help)
      echo "Arguments to give:"
      echo "-gff  :         path of the annotation file"
      echo "-fasta:         path of the genome file"
      echo "-fastq:         path of the riboseq file"
      echo "-adapt:         adaptor sequence"
      echo "-kmer :         range of reads size to analyse"
      exit 0
      ;;
  esac
done;
#Launch python scripts
python ORFphase_python3.py -gff "${GFF}" -fasta "${FASTA}" -fastq "${FASTQ}" -adapt "${ADAPT}" -kmer kmer1 kmer2
j=$(python find_best_kmer.py -fastq "${FASTQ}")
python ORFribomap.py -gff "${GFF}" -fasta "${FASTA}" -fastq "${FASTQ}" -kmer "$j" -cutdir . -o nc -features_include nc_intergenic
