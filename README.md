# TCG
Here are the codes related to "A core microbiome signature as an indicator of health" paper (DOI: 10.1016/j.cell.2024.09.019)

#1.	Data quality control of metagenomes
Tool: kneaddata (https://huttenhower.sph.harvard.edu/kneaddata/)
Input: $R1: R1 reads, $R2: R2 reads
kneaddata -i1 $R1 -i2 $R2 -o $prefix -db $refdb --output-prefix $prefix \
-t 20 --decontaminate-pairs strict --sequencer-source $sequencer_source --run-trim-repetitive \
--trimmomatic-options=\"ILLUMINACLIP:$adapter:2:30:10 SLIDINGWINDOW:4:20 MINLEN:60\" \
--run-fastqc-start --run-fastqc-end --bypass-trf \
--trimmomatic Trimmomatic-0.39 \\
--fastqc FastQC --remove-intermediate-output –reorder
