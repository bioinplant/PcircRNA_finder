# PcircRNA_finder
PcircRNA_finder: a software for circRNA prediction in plants

Recent studies in animals have found a kind of non-coding RNA molecule called circular RNA (circRNA) produced by splicing complex with one or more exons of the gene. In animals and human cells, circular RNA has miRNA sponge function, regulating its parent gene functions and also it was reported that the formation of circular RNA depends on both sides of the intron complementary sequence. However, the situation of circular RNA in plants is unknown. Using the current softwares or algorithms developed for animals or humans, the accuracy of the prediction of circular RNA in plant is relatively low and the sensitivity is not high. Here, a new circRNA prediction method called PcircRNA_finder is developed in plants. Using the simulated data and real RNA-seq data, we show that the software has good sensitivity, and can be more comprehensive to find out the circular RNA. It can be used to distinguish and predict circular RNA with alternative donor/acceptor sites with minor changes. At the same time, the software has a better output, and gives a more detailed information on the circular RNA, which facilitate the follow-up of experimental verification and the functional study of circular RNA.

 

	PcircRNA_finder is mainly designed for exonic circRNA prediction (Figure 1), which consists of three parts: 
	(1) Catchor: Catchor is used to collect as many backsplice sites as possible, which are supported by chiastic clipping mapping based on many available fusion detection methods, including Tophat-Fusion (Kim and Salzberg, 2011), STAR-Fusion (Dobin, et al., 2013), find_circ (Memczak et al., 2013), Mapsplice (Wang, et al., 2010), segemehl (Hoffmann, et al., 2014). Among these candidate back-splice sites, it has many false positive sites for specific sofware and then, depends on next Filter module to filter it out.
 	(2) Annotator: Annotator can be used to annotate the candidate sites with gene structure annotation, which compares the back-splice sites to the exon-intron boundary. And It's based on two experimental observations, one is that the experimental verified circRNAs (10/18) in rice are consistent with exon boundary (Ye, et al., 2015), which shows circRNA used same splicing site with pre-mRNA splicing; the other is that the biogensis of circRNA is dependent on canonical spliceosomal machinery even though circRNA's back-splicing site is cryptic and flexible (Starke, et al., 2015). That's to say, circRNA is another form of alternative splicing. 
	(3) Filter: Filter is used to prioritize the above candidate circRNA. It has two parts: filtering the candidate circRNAs with two kinds of splicing signal, namely U2 based spliceosome (usually with consensus sequence GT-AG and GC-AG ) and U12 based minor spliceosome (usually with consensus sequence AT-AC) (Reddy, et al., 2013; Staiger and Brown, 2013); creating a pseudoRef with chiastic backsplice site flanking sequences and then mapping raw reads to it to check the back-splice sites.

Installation Prerequisite
1. Python (version=2.7)
https://www.python.org/
2.perl (version=5.10)
https://www.perl.org/
3. STAR
https://github.com/alexdobin/STAR
4. tophat
http://ccb.jhu.edu/software/tophat/fusion_index.shtml
5. find_circ
http://circbase.org/download/find_circ.tar.gz
6. mapsplice
http://www.netlab.uky.edu/p/bioinfo/MapSplice2
7. segemehl
http://www.bioinf.uni-leipzig.de/Software/segemehl/
8. PcircRNA_finder
9. bedtools
https://sourceforge.net/projects/bedtools/
10. samtools
http://samtools.sourceforge.net/
11. bowtie2
http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
12. bowtie1
http://bowtie-bio.sourceforge.net/index.shtml

      Running step
      1. Running STAR to get "stiChimeric.out.junction" file, example command below:
      STAR --runMode genomeGenerate --genomeDir /home3/cl/at_sti_15x/STAR_index --genomeFastaFiles Athaliana_167.fa  --runThreadN 8

      STAR --genomeDir  /home3/cl/at_sti_15x/STAR_index  --readFilesIn ../left.fastq  ../right.fastq --runThreadN 8 --sjdbGTFfile /home1/cl/ok/TAIR10_GFF3_genes.gtf  --chimSegmentMin 20 --chimScoreMin 1 --alignIntronMax 100000 --outFilterMismatchNmax 4 --alignTranscriptsPerReadNmax 100000 --outFilterMultimapNmax 2 --outFileNamePrefix   sti  --outSAMtype BAM SortedByCoordinate

      2. Running tophat to get " accepted_hits.bam " file, example command below:
      tophat2 -a 6 --microexon-search -m 2 -p 10 -G  $gtf_file   -o ${sample_name}   $genome_bowtie_index $fastq1list  $fastq2list

      bamToFastq -i   ${sample_name}/unmapped.bam  -fq  ${sample_name}/unmapped.fastq

      tophat2 -o   ${sample_name}_fusion -p 15 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search   $genome_bowtie_index      ${sample_name}/unmapped.fastq

      3. Running find_circ to get " find_circ_circRNA.bed" file, example command below:
      $bowtie2_path  -p16 --very-sensitive --mm -M20 --score-min=C,-15,0 -x  $genome_index  -q  -U  $fastq_file   2> ${name}_bt2_firstpass.log | $samtools_path  view -hbuS - | $samtools_path  sort -  $name

      $samtools_path   view -hf 4  ${name}.bam | $samtools_path  view -Sb - > unmapped_${name}.bam

      python  $find_circRNA_path/unmapped2anchors.py   unmapped_${name}.bam  > unmapped_${name}_anchors.qfa

      $bowtie2_path  --reorder --mm -M20 --score-min=C,-15,0 -q -x  $genome_index  -U unmapped_${name}_anchors.qfa  2> bt2_secondpass_test.log | python  $find_circRNA_path/find_circ.py -r $sample_description_file  -G $chr_dir -p $prefix -s test_out_$name/sites.log > test_out_$name/sites.bed 2> test_out_$name/sites.reads

      grep  _circ_   sites.bed    > find_circ_circRNA.bed

      4. Running mapsplice to get " fusions_raw.txt " file, example command below:
      python  /datacenter/disk3/cl/mapsplice/MapSplice-v2.1.8/mapsplice.py   -p 10  --non-canonical  --fusion-non-canonical   --min-fusion-distance  200  -c  /datacenter/disk1/cl/data/at/chr    -x  /home1/cl/ok/Athaliana_167.fa   --gene-gtf  /home1/cl/ok/TAIR10_GFF3_genes.gtf     -1  ../left.fastq      -2    ../right.fastq    --qual-scale  phred33   -o      circ  

      5. Running segemehl to get " splicesites.bed" file, example command below:
      $segemehl_path/segemehl.x -t 8 -s -d $genome_fasta -i $genome_segemehl_index  -q $fastq1_file   -p $fastq2_file -o ${name_1}.sam -u ${name_1}.fq -S -T -D 2

      samtools view -hbuS   ${name_1}.sam | samtools sort - ${name_1}.sorted

      samtools view -h ${name_1}.sorted.bam | gzip -c > ${name_1}.sorted.sam.gz

      $segemehl_path/testrealign.x -d $genome_fasta  -q   ${name_1}.sorted.sam.gz -U ${name_1}.splitmap.txt -T ${name_1}.transmap.txt -n

      6. Running PcircRNA_finder to get final results, example command below:
      perl ecircRNA_finder.pl /home3/cl/at_sti_15x/re_run1/1/1/stiChimeric.out.junction  /home3/cl/at_sti_15x/re_run1/1/1/sti_fusion/accepted_hits.bam /home3/cl/at_sti_15x/re_run1/1/1/test_out_sti/find_circ_circRNA.bed  /home3/cl/at_sti_15x/re_run1/1/1/circ/fusions_raw.txt  /home3/cl/at_sti_15x/re_run1/1/1/splicesites.bed   20000  /home3/cl/at_sti_15x/re_run1/1/1/at.txt 5 5  /home1/cl/ok/Athaliana_167.fa  100  1  left.fastq right.fastq  1


Output files
column 1: canonical splicing corresponding position
column 2: circRNA position
column 3: chromosome 
column 4: start (0-based)
column 5: end
column 6: name
column 7: -
column 8: strand
column 9:host gene &circRNA's included exons
column 10: Canonical means canonical splicing sites are the same as circRNA position; Alt_acceptor means circRNA position have alternative acceptor sites compared to canonical splicing sites; Alt_donor means circRNA position have alternative donor sites compared to canonical splicing sites
column 11: backspliced reads number

example output file: left_Final_exonic_circRNA.txt
Chr2:9845487..9848769	Chr2:9845486..9848769	Chr2	9845486	9848769	backjunc_117	0,16,30,17,0	-	AT2G23140.1:exon_5,exon_4,exon_3,exon_2||AT2G23140.2:exon_5,exon_4,exon_3,exon_2	Alt_acceptor	16
Chr3:10084136..10086696	Chr3:10084136..10086696	Chr3	10084135	10086696	backjunc_118	0,19,0,23,0	-	AT3G27300.1:exon_12,exon_11,exon_10,exon_9,exon_8,exon_7,exon_6,exon_5,exon_4,exon_3,exon_2,exon_1||AT3G27300.2:exon_13,exon_12,exon_11,exon_10,exon_9,exon_8,exon_7,exon_6,exon_5,exon_4,exon_3,exon_2,exon_1||AT3G27300.3:exon_12,exon_11,exon_10,exon_9,exon_8,exon_7,exon_6,exon_5,exon_4,exon_3,exon_2,exon_1	Canonical	6
Chr3:10084137..10086697	Chr3:10084136..10086696	Chr3	10084136	10086697	backjunc_119	0,0,36,0,0	-	AT3G27300.1:exon_12,exon_11,exon_10,exon_9,exon_8,exon_7,exon_6,exon_5,exon_4,exon_3,exon_2,exon_1||AT3G27300.2:exon_13,exon_12,exon_11,exon_10,exon_9,exon_8,exon_7,exon_6,exon_5,exon_4,exon_3,exon_2,exon_1||AT3G27300.3:exon_12,exon_11,exon_10,exon_9,exon_8,exon_7,exon_6,exon_5,exon_4,exon_3,exon_2,exon_1	Alt_acceptor,Alt_donor	12
