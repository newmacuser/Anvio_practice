#!/bin/bash

echo -e "Please enter path to working dir" 
read wdir 
cd $wdir 
echo $PWD 
echo -e "Please enter path to raw sequencing files" 
read seqdir
echo $seqdir 
echo -e "Please enter sample name" 
read name 
echo -e "Please make sure you have the following files in working dir: 'sample_name.contig.fasta' 'sample.txt' 'sample_list.txt' 'contigs_split.R' 'bin_naming_wsl_refinem.R'" 
echo -e "Please enter yes to continue"
read yes
if [[ $yes = yes ]]
then
  echo "Script continue"
else
  echo "Stop and check your files"
  exit 1
fi
echo -e "Please specify numbers of threads" 
read thread
echo -e "Please specify full path (include file name) to kaiju database 'nodes.dmp' file" 
read knode
echo -e "Please specify full path (include file name) to kaiju database 'names.dmp' file" 
read kname
echo -e "Please specify full path (include file name) to kaiju database 'kaiju_db_refseq.fmi' file" 
read kfmi
echo -e "Please specify full path (include file name) to refinem database 'gtdb_r89_protein_db.2019-09-27.faa.dmnd' file" 
read dmnd
echo -e "Please specify full path (include file name) to refinem database 'gtdb_r89_taxonomy.2019-09-27.tsv' file" 
read taxo
echo -e `date '+%Y-%m-%d %H:%M:%S %A'` Script starting



mkdir 03_CONTIGS
anvi-script-reformat-fasta $name.contigs.fasta -o 03_CONTIGS/contigs.fa --min-len 2500 --simplify-names --report name_conversions.txt
mkdir 04_MAPPING
bowtie2-build 03_CONTIGS/contigs.fa 04_MAPPING/contigs

for sample in `awk '{print $1}' sample.txt`
do
if [ "$sample" == "sample" ]; then continue; fi
R1s=`ls $seqdir/$sample*_1.fastq | python2 -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
R2s=`ls $seqdir/$sample*_2.fastq | python2 -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
echo $sample
echo $R1s
echo $R2s
bowtie2 --threads $thread -x 04_MAPPING/contigs -1 $R1s -2 $R2s --no-unal -S 04_MAPPING/$sample.sam
samtools view --threads $thread -F 4 -bS 04_MAPPING/$sample.sam > 04_MAPPING/$sample-RAW.bam
anvi-init-bam 04_MAPPING/$sample-RAW.bam -o 04_MAPPING/$sample.bam
rm 04_MAPPING/$sample.sam 04_MAPPING/$sample-RAW.bam
done

rm 04_MAPPING/*.bt2

mkdir 05_contigsDatabase
anvi-gen-contigs-database -f 03_CONTIGS/contigs.fa -o 05_contigsDatabase/contigs.db -n $name 

for i in $(cat  sample_list.txt)
 do
echo $i
anvi-profile -i 04_MAPPING/"$i".bam -c 05_contigsDatabase/contigs.db -T $thread --min-contig-length 2500 --output-dir 06_profiling/"$i" --sample-name "$i"
echo Finish "$i"
 done

cd 06_profiling
anvi-merge */PROFILE.db -o merged_profile -c ../05_contigsDatabase/contigs.db --enforce-hierarchical-clustering
cd ..
anvi-export-splits-and-coverages -p 06_profiling/merged_profile/PROFILE.db -c 05_contigsDatabase/contigs.db -o export_splites -O $name --report-contigs
Rscript contigs_split.R
echo -e `date '+%Y-%m-%d %H:%M:%S %A'` Finish mapping and profiling

# Run cogs and kaiju
echo -e `date '+%Y-%m-%d %H:%M:%S %A'` Start to run COGs and Kaiju
anvi-run-hmms -c 05_contigsDatabase/contigs.db --num-threads $thread
anvi-run-ncbi-cogs -c 05_contigsDatabase/contigs.db  -T $thread
anvi-get-sequences-for-gene-calls -c 05_contigsDatabase/contigs.db -o 05_contigsDatabase/gene_calls.fa
cd 05_contigsDatabase
kaiju -t $knode -f $kfmi -i gene_calls.fa -o gene_calls_ref.out -z $thread -v 
kaiju-addTaxonNames -t $knode -n $kname -i gene_calls_ref.out -o gene_calls.names -r superkingdom,phylum,order,class,family,genus,species 
anvi-import-taxonomy-for-genes -c contigs.db -i gene_calls.names -p kaiju --just-do-it 
cd ..
echo -e `date '+%Y-%m-%d %H:%M:%S %A'` Finish COGs and Kaiju

# Run Maxbin2 and Refinem
echo -e `date '+%Y-%m-%d %H:%M:%S %A'` Start to run Maxbin2 and Refinem
mkdir 07_maxbin
cd 07_maxbin
run_MaxBin.pl -contig ../export_splites/$name-CONTIGS.fa -out maxbin2 -thread $thread -abund ../export_splites/abund.1.txt -abund2 ../export_splites/abund.2.txt -abund3 ../export_splites/abund.3.txt -abund4 ../export_splites/abund.4.txt -abund5 ../export_splites/abund.5.txt
mkdir refinem 
cd refinem 
refinem scaffold_stats -c $thread -x fasta ../../export_splites/R1-CONTIGS.fa ../ ./stats ../../04_MAPPING/$name.bam 
refinem outliers ./stats/scaffold_stats.tsv ./outliers 
refinem filter_bins -x fasta ../  ./outliers/outliers.tsv ./filtered_bins 
refinem call_genes -c $thread -x fasta ./filtered_bins ./gene_output 
refinem taxon_profile -c $thread ./gene_output ./stats/scaffold_stats.tsv $dmnd $taxo ./taxon_profile 
refinem taxon_filter -c $thread  ./taxon_profile taxon_filter.tsv 
refinem filter_bins -x fasta ./filtered_bins taxon_filter.tsv ./tax_filtered_bins 
cd tax_filtered_bins 
for i in *.filtered.filtered.fasta 
do 
grep ">" $i > ${i%.*}.txt 
done 
Rscript ../../../bin_naming_wsl_refinem.R
cd ../../.. 
anvi-import-collection  07_maxbin/refinem/tax_filtered_bins/bins.txt -p 06_profiling/merged_profile/PROFILE.db -c 05_contigsDatabase/contigs.db -C maxbin --contigs-mode 
echo -e `date '+%Y-%m-%d %H:%M:%S %A'` Finish Maxbin2 and Refinem
anvi-summarize -p 06_profiling/merged_profile/PROFILE.db  --report-aa-seqs-for-gene-calls -c 05_contigsDatabase/contigs.db -o 08_summary -C maxbin 
echo -e `date '+%Y-%m-%d %H:%M:%S %A'` All finished!
