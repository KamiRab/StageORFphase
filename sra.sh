#!/usr/bin/env zsh
set -euxo pipefail
arr_srr=($(tail -n +1 PolyAribo.csv | cut -d ',' -f2 ))
for i in "${arr_srr[@]}"
do
  if [[ -f "$i".fastq ]];then
    echo "$i.fastq on external drive"
  elif [ -f /media/crabier/Elements/SRA/"$i".fastq ];then
      echo "$i.fastq exists already on ext"
  else
    # prefetch "$i"
    fasterq-dump "$i" || continue
    mv "$i".fastq 
  fi
done
arr_srr2=($(tail -n +1 TrueAdapt.csv | cut -d ',' -f2 ))
for i in "${arr_srr2[@]}"
do
  if [ -f /media/crabier/Elements/SRA/"$i".fastq ];then
    echo "$i.fastq exists already"
  else
    prefetch "$i"
    fasterq-dump "$i" -O /media/crabier/Elements/SRA || continue
  fi
done
# sudo wondershaper wlp1s0 256 128
# for i in "${arr_srr[@]}"
# do
#   scp /media/crabier/Elements/SRA/"$i".fastq i2bc:/store/EQUIPES/BIM/MEMBERS/camille.rabier
# done
# for i in "${arr_srr2[@]}"
# do
#   scp /media/crabier/Elements/SRA/"$i".fastq i2bc:/store/EQUIPES/BIM/MEMBERS/camille.rabier
# done
