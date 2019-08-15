# The UPHL-Reference-Free pipeline takes paired-end fastq files to contigs for microbial WGS.

Below is a list of the commands used in UPHL's Reference-Free Workflow for those who wish to use the same commands as UPHL, but with their own workflow manager. There are also custom scripts to format the results for [MultiQC](https://github.com/ewels/MultiQC) custom editions. As time permits, we hope to add to MultiQC's supported tools and remove our personal custom scripts. 

- [seqyclean](https://github.com/ibest/seqyclean)  
```
seqyclean -minlen 25 -qual -c /Adapters_plus_PhiX_174.fasta -1 Sequencing_reads/Raw/sample_1.fastq -2 Sequencing_reads/Raw/sample_2.fastq -o Sequencing_reads/QCed/sample_clean
```
- [shovill](https://github.com/tseemann/shovill)
```
shovill --cpu 1 --ram $RAM --outdir shovill_result/sample --R1 Sequencing_reads/QCed/sample_clean_PE1.fastq --R2 Sequencing_reads/QCed/sample_clean_PE2.fastq
```
- [prokka](https://github.com/tseemann/prokka)
```
prokka --cpu 1 --compliant --centre --URF --mincontiglen 500 --outdir Prokka/sample --locustag locus_tag --prefix sample --genus ${mash_result[0]} --species ${mash_result[1]} --force shovill_result/sample/contigs.fa
```
- [fastqc](https://github.com/s-andrews/FastQC)
```
fastqc --outdir fastqc --threads 1 Sequencing_reads/*/*.fastq*
```
- [cg-pipeline](https://github.com/lskatz/CG-Pipeline)
```
run_assembly_shuffleReads.pl -gz Sequencing_reads/QCed/sample_clean_PE1.fastq Sequencing_reads/QCed/sample_clean_PE2.fastq > Sequencing_reads/shuffled/sample_clean_shuffled.fastq.gz
run_assembly_readMetrics.pl Sequencing_reads/shuffled/sample_clean_shuffled.fastq.gz --fast --numcpus 1 -e $genome_length
```
- [quast](https://github.com/ablab/quast)
```
quast.py ALL_assembled/sample_contigs.fa --output-dir quast/sample --threads 1
```
- [multiqc](https://github.com/ewels/MultiQC)
```
multiqc -f --outdir logs --cl_config "prokka_fn_snames: True" .
```
- [mash](https://github.com/marbl/Mash)
```
cat Sequencing_reads/QCed/sample_clean_PE1.fastq Sequencing_reads/QCed/sample_clean_PE2.fastq | mash sketch -m 2 -o mash/sample -
mash dist -p 1 -v 0 /db/RefSeqSketchesDefaults.msh mash/sample.msh | sort -gk3 > mash/sample_mashdist.txt
```
- [seqsero](https://github.com/denglab/SeqSero)
```
SeqSero.py -m 2 -d SeqSero/sample -i Sequencing_reads/QCed/sample_clean_PE1.fastq Sequencing_reads/QCed/sample_clean_PE2.fastq
```
- [abricate](https://github.com/tseemann/abricate)
```
abricate --db serotypefinder --threads 1 shovill_result/sample/contigs.fa > abricate_results/serotypefinder/serotypefinder.sample.out.tab
abricate --summary abricate_results*/serotypefinder/serotypefinder*tab
abricate --db ncbi --threads 1 shovill_result/sample/contigs.fa > abricate_results/ncbi/ncbi.sample.out.tab
abricate --summary abricate_results*/ncbi/ncbi*tab
abricate --db vfdb --threads 1 shovill_result/sample/contigs.fa > abricate_results/vfdb/vfdb.sample.out.tab
abricate --summary abricate_results*/vfdb/vfdb*tab > abricate_results/vfdb/vfdb.summary.txt
```
- [blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
```
blastn -query shovill_result/sample/contigs.fa -out blast/sample.tsv -num_threads 1 -db /blast/blastdb/nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 10 -max_hsps 1 -evalue 1e-25
```
- [bwa](http://bio-bwa.sourceforge.net/)
```
bwa index shovill_result/sample/contigs.fa
bwa mem -t 1 shovill_result/sample/contigs.fa Sequencing_reads/QCed/sample_clean_PE1.fastq Sequencing_reads/QCed/sample_clean_PE2.fastq | samtools sort -o bwa/sample.sorted.bam
```
- [blobtools](https://blobtools.readme.io/docs)
```
blobtools create -o blobtools/sample -i shovill_result/sample/contigs.fa -b bwa/sample.sorted.bam -t blast/sample.tsv
blobtools view -i blobtools/sample.blobDB.json -o blobtools/
blobtools plot -i blobtools/sample.blobDB.json -o blobtools/ -r species --format png
```
