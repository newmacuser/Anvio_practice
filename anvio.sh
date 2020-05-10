#!/bin/bash

echo -e "Please enter path of working dir" 
read wdir 
cd $wdir 
echo $PWD 
echo -e "Please enter sample name" 
read name 
echo -e "Please make sure you have the following files in working dir: 'sample_name.contig.fasta' 'sample.txt' 'sample_list.txt'" 
echo -e "Please enter yes to continue"
read yes
if [[ $yes = yes ]]
then
  echo "Script starting"
  echo -e `date '+%Y-%m-%d %H:%M:%S %A'`
else
  echo "Stop and check your files"
  exit 1
fi

mkdir 03_CONTIGS
anvi-script-reformat-fasta $name.contigs.fasta -o 03_CONTIGS/contigs.fa --min-len 2500 --simplify-names --report name_conversions.txt
mkdir 04_MAPPING
bowtie2-build 03_CONTIGS/contigs.fa 04_MAPPING/contigs

for sample in `awk '{print $1}' sample.txt`
do
if [ "$sample" == "sample" ]; then continue; fi
R1s=`ls /mnt/e/RWSAX/metagenomic_rawdata/$sample*_1.fastq | python2 -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
R2s=`ls /mnt/e/RWSAX/metagenomic_rawdata/$sample*_2.fastq | python2 -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
echo $sample
echo $R1s
echo $R2s
bowtie2 --threads 20 -x 04_MAPPING/contigs -1 $R1s -2 $R2s --no-unal -S 04_MAPPING/$sample.sam
samtools view --threads 12 -F 4 -bS 04_MAPPING/$sample.sam > 04_MAPPING/$sample-RAW.bam
anvi-init-bam 04_MAPPING/$sample-RAW.bam -o 04_MAPPING/$sample.bam
rm 04_MAPPING/$sample.sam 04_MAPPING/$sample-RAW.bam
done

rm 04_MAPPING/*.bt2

mkdir 05_contigsDatabase
anvi-gen-contigs-database -f 03_CONTIGS/contigs.fa -o 05_contigsDatabase/contigs.db -n $name 
anvi-run-hmms -c 05_contigsDatabase/contigs.db --num-threads 12
anvi-run-ncbi-cogs -c 05_contigsDatabase/contigs.db  -T 12
anvi-get-sequences-for-gene-calls -c 05_contigsDatabase/contigs.db -o 05_contigsDatabase/gene_calls.fa
cd 05_contigsDatabase
centrifuge -f -x /mnt/e/Database/centrifuge_db/p_compressed gene_calls.fa -S centrifuge_hits.tsv -p 12
wc -l centrifuge_report.tsv centrifuge_hits.tsv
anvi-import-taxonomy-for-genes -c contigs.db -i centrifuge_report.tsv centrifuge_hits.tsv -p centrifuge

cd ..
for i in $(cat  sample_list.txt)
 do
echo $i
anvi-profile -i 04_MAPPING/"$i".bam -c 05_contigsDatabase/contigs.db -T 16 --min-contig-length 2500 --output-dir 06_profiling/"$i" --sample-name "$i"
echo Finish "$i"
 done

cd 06_profiling
anvi-merge */PROFILE.db -o merged_profile -c ../05_contigsDatabase/contigs.db --enforce-hierarchical-clustering
cd ..
anvi-export-splits-and-coverages -p 06_profiling/merged_profile/PROFILE.db -c 05_contigsDatabase/contigs.db -o export_splites -O $name --report-contigs
echo -e `date '+%Y-%m-%d %H:%M:%S %A'` Finish
