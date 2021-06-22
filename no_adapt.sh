#!/usr/bin/env zsh
set -euxo pipefail
arr_srr=(SRR1520311 SRR1520312 SRR1520313 SRR1520314 SRR1520315 SRR1520316 SRR1520317 SRR1520318 SRR1520319 SRR1520320 SRR1520321 SRR1520322 SRR1520323 SRR1520324 SRR1520325 SRR1520326 SRR1520327 SRR1520328 SRR1520329 SRR1520330 SRR1520331 SRR1520332 SRR1520333 SRR1520334
)
for i in "${arr_srr[@]}"
do
  if [[ -f ./SRR_to_treat/"$i".fastq ]];then
    for j in {26..30}
      cutadapt ./SRR_to_treat/"$i".fastq -j 1 --quality-base=33 -m "$j" -M "$j" -e 0.12 -o ./kmer_"$j"/"$i"_kmer_"$j".fastq
  else
    # prefetch "$i"
    echo "No fastq"
  fi
done
