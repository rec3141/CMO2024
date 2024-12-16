#!/bin/bash
# quick look at nanopore data

cd /data/project_CMO

# quality control with bbduk
for file in CMO-20241*.fastq.gz; do ~/apps/bbmap/bbduk.sh in=$file out=${file}_qc.fq.gz minlength=1000 minavgquality=20; done

# do flye assembly at qc20
cat CMO*_qc.fq.gz > CMO-all.qc.fq.gz
~/apps/bbmap/reformat.sh in=CMO-all.qc.fq.gz out=CMO-all.qc.fa.gz

flye --nano-raw CMO-all.qc.fq.gz --out-dir /data/project_CMO/flye_qc20_dec14 --threads 12 --meta --min-overlap 1000

# in Bandage: save entire graph to FASTA: all_positive_graph_nodes.fasta

#qc15 was too big for Grid to fit in 100Gb memory, 30Gb of seqs
#qc18 probably too big, almost filled memory, 10Gb of seqs
#qc20 was 1Gb of seqs

# do another flye assembly at qc18
for file in CMO-20241*.fastq.gz; do ~/apps/bbmap/bbduk.sh in=$file out=${file}_qc18.fq.gz minlength=1500 minavgquality=18; done

cat CMO*_qc18.fq.gz > CMO-all.qc18.fq.gz
~/apps/bbmap/reformat.sh in=CMO-all.qc18.fq.gz out=CMO-all.qc18.fa.gz

flye --nano-raw CMO-all.qc18.fq.gz --out-dir /data/project_CMO/flye_qc18_dec14 --threads 12 --meta --min-overlap 1000

#flye --nano-raw CMO-all.qc18.fq.gz --out-dir /data/project_CMO/flye_qc18_2k_dec14 --threads 12 --meta --min-overlap 2000

#flye --nano-raw CMO-all.qc18.fq.gz --out-dir /data/project_CMO/flye_qc18_5k_dec14 --threads 12 --meta --min-overlap 5000

# sketch all reads with fairy
for file in CMO*fastq.gz; do fairy sketch -r $file -d fairy_all_out; done

# sample all coverage with fairy
fairy coverage fairy_qc20_out/*.bcsp /data/project_CMO/flye_qc20_dec14/all_positive_graph_nodes.fasta -t 12 -o fairy_qc20_out/assembly_graph_noqc_cov.tsv


# sketch QC reads with fairy
for file in CMO*_qc.fq.gz; do fairy sketch -r $file -d fairy_qc20_out; done
for file in CMO*_qc18.fq.gz; do fairy sketch -r $file -d fairy_qc18_out; done

# sample QC coverage with fairy
fairy coverage fairy_qc20_out/*.bcsp /data/project_CMO/flye_qc20_dec14/all_positive_graph_nodes.fasta -t 12 -o fairy_qc20_out/assembly_graph_cov.tsv
fairy coverage fairy_qc18_out/*.bcsp /data/project_CMO/flye_qc18_dec14/all_positive_graph_nodes.fasta -t 12 -o fairy_qc18_out/assembly_graph_cov.tsv

~/apps/metabat/bin/metabat2 -i /data/project_CMO/flye_qc20_dec14/all_positive_graph_nodes.fasta -o /data/scratch/fairy_qc20_out/metabat_graph_out -m 1500 -s 60000 --saveCls --unbinned -v -a /data/scratch/fairy_qc20_out/assembly_graph_cov.tsv 
#[00:00:01] 82.54% (57506759 bases) of large (>=1500) and 11.17% (25837 bases) of small (<1500) contigs were binned.
# 139 bins (57532596 bases in total) formed.

~/apps/metabat/bin/metabat2 -i /data/project_CMO/flye_qc18_dec14/all_positive_graph_nodes.fasta -o /data/scratch/fairy_qc18_out/metabat_graph_out -m 1500 -s 0 --saveCls --unbinned -v -a /data/scratch/fairy_qc18_out/assembly_graph_cov.tsv 


cd fairy_qc20_out
cd fairy_qc18_out

# taxonomic ids using sendsketch
for file in `ls -S metabat_graph_out*.fa`; do ~/apps/bbmap/sendsketch.sh in=$file address=protein translate=t out=$file.sk.tsv format=2 printall name0=$file printname0=t ow; done

# TNF calculation
# run on bins
grep '>' /data/project_CMO/flye_qc20_dec14/all_positive_graph_nodes.fasta | sed 's/>//' | paste - - -  > graph.annotation.txt
perl ~/apps/tetramer_freqs_esom.pl -f /data/project_CMO/flye_qc20_dec14/all_positive_graph_nodes.fasta -a graph.annotation.txt -min 1500 -max 10000000
head -n4 Tetra_all_positive_graph_nodes_1500.lrn | tail -n1 | cut -f2- -d' ' > tnfs.txt

grep '>' /data/project_CMO/flye_qc18_dec14/all_positive_graph_nodes.fasta | sed 's/>//' | paste - - -  > graph.annotation.txt
perl ~/apps/tetramer_freqs_esom.pl -f /data/project_CMO/flye_qc18_dec14/all_positive_graph_nodes.fasta -a graph.annotation.txt -min 1500 -max 10000000
head -n4 Tetra_all_positive_graph_nodes_1500.lrn | tail -n1 | cut -f2- -d' ' > tnfs.txt

# run on assembly, split into 5k bp chunks 
cat /data/scratch/fairy_qc18_out/metabat_graph_out*.fa > metabat_graph_out_all.fasta
grep '>' metabat_graph_out_all.fasta | sed 's/>//' | paste - - -  > annotation.graph_1500.txt
perl ~/apps/tetramer_freqs_esom.pl -f metabat_graph_out_all.fasta -a annotation.graph_1500.txt -min 1500 -max 5000

# bin split coverage with fairy QC
fairy coverage /data/scratch/fairy_qc18_out/*qc18*.bcsp Tetra_metabat_graph_out_all_1500_5000_split.fasta -t 12 -o metabat_graph_cov.tsv

# bin split coverage with fairy noQC
fairy coverage $(find /data/scratch/ -name "*.bcsp" | grep -v qc) Tetra_metabat_graph_out_all_1500_5000_split.fasta -t 12 -o metabat_graph_noqc_cov.tsv

# run TNF on reads
zgrep '>' ./../CMO-all.qc.fa.gz | sed 's/>//' | paste - - -  > reads.annotation.txt
perl ~/apps/tetramer_freqs_esom.pl -f <(gunzip -c ./../CMO-all.qc.fa.gz) -a reads.annotation.txt -min 1500 -max 10000000


# reads coverage with fairy
fairy coverage /data/scratch/fairy_qc18_out/*qc18*.bcsp /data/scratch/fairy_qc20_out/Tetra_63_1500_10000000_split.fasta -t 12 -o reads_cov.tsv

# add filenames to bin contigs
for file in metabat_graph_out.*.fa; do grep '>' $file | cut -f2 -d'>' | cut -f1 -d$'\t' | awk -v fname=$file '{ print fname "\t" $0 }'; done  > metabat_binids.tsv

for file in metabat_graph_out.*.fa.sk.tsv; do grep -H '%' -m1 $file; done > besthits.sk.tsv
cut -f1,37 -d$'\t' besthits.sk.tsv > tophits.sk.tsv


