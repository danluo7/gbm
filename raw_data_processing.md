# gbm

## check quality of reads

    fastqc *.fastq.gz 
    
    python3 -m multiqc . 
    
    mkdir fastqc
    mv *fastqc* fastqc
    
hg38 genome index:    
  
Sources for obtaining gene annotation files formatted for HISAT2: HISAT2 pre-built Genome Index 
http://daehwankimlab.github.io/hisat2/download/

Also available from their FTP site ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/
(this was used by Andrew in Hingtgen's lab)

    wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/hg38.tar.gz
    tar -xzvf hg38.tar.gz
  
These 8 files together constitute the index: they are all that is needed to align reads to that reference. 


For mouse mm10 genome hisat2 pre-built index:

    mkdir -p $gbm/RNA_REF_FA
    cd $gbm/RNA_REF_FA
    wget https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz
    tar -xzvf mm10_genome.tar.gz

note: GRCm39 = mm39, and GRCm38 = mm10



Also download the reference genome compatible to Hisat2 from UCSC: 
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/
for other genomes: https://hgdownload.soe.ucsc.edu/downloads.html

File used by Andrew is:

    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz 
    gzip -d hg38.ncbiRefSeq.gtf.gz


For mouse reference UCSC mm10: https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/

    mkdir -p $gbm/RNA_REF_GTF
    cd $gbm/RNA_REF_GTF
    wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ncbiRefSeq.gtf.gz
    gzip -d mm10.ncbiRefSeq.gtf.gz

note: mouse reference UCSC mm39 became available in 2020 (mm38/mm10 came out in 2011) but not used by hisat2: https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/
reason this matters: mm10 and GRCm38 are synonymous: https://github.com/kundajelab/atac_dnase_pipelines/issues/143 https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.20/ except: UCSC version will have chr* identifiers in the row names. 


Alternatively, UCSC's references are also hosted on S3 bucket: http://daehwankimlab.github.io/hisat2/download/



### alignment using hisat2 (which suceeded hisat and tophat)

setting the Read Group info:

a) Sample and library tags. Can be autmatically detected from current sample naming scheme:
<Sample.ID><Index.Sequence><Lane.ID><Set.number>.fastq

which for this dataset is: 
6931_1_S1_L001_R1_001.fastq.gz

sample.ID: 6931
index.sequence: 1
lane.id: L001
set.number: R1


Also can be identified from the name of a sequence read in the Fastq file:

@(instrument id) @A00201R
:(run number) 402
:(flowcell ID) HCTMMDRXY
:(lane) 1
:(tile) 2102
:(x_pos) 21305
:(y_pos) (read) 28087
:(pair) 2
:(is filtered) N
:(control number) 0
:(UMI) GCTAATAGGA+AGAGCACTAG
:(index sequence)

ID and PU (to enable merging replictes)


ID = Read group identifier = {FLOWCELL BARCODE/ID}.{LANE} HCTMMDRXY.1 or .2
SM = <Sample.ID> 1 through xx (individual sample that was submitted)
LB = <Sample.ID>_<Index.Sequence> (to indentify if same library was used)

PU = Platform Unit = {FLOWCELL_BARCODE}.{LANE}.{library-specific identifier}. This is the most specific definition for a group of reads.

## At this point decide on naming conventions. Order the samples logically 

(ie. if samples 1 2 3 should be ordered 3 1 2, then in hisat2 put input sample file 1 > 3.bam, and sample file 2 > 1.bam etc. Then make a note of this change.


	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:1 --rg LB:1_TGGTAGAGAT+TGTTGTTCGT --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.TGGTAGAGAT+TGTTGTTCGT -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_1_S1_L001_R1_001.fastq.gz -2 $gbm_data/6931_1_S1_L001_R2_001.fastq.gz -S $gbm_data/alignments/1_rep1.sam
	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:1 --rg LB:1_TGGTAGAGAT+TGTTGTTCGT --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.TGGTAGAGAT+TGTTGTTCGT -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_1_S1_L002_R1_001.fastq.gz -2 $gbm_data/6931_1_S1_L002_R2_001.fastq.gz -S $gbm_data/alignments/1_rep2.sam

	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:2 --rg LB:2_AGTACTCATG+GTAGAGTCAG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.AGTACTCATG+GTAGAGTCAG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_2_S2_L001_R1_001.fastq.gz -2 $gbm_data/6931_2_S2_L001_R2_001.fastq.gz -S $gbm_data/alignments/2_rep1.sam
	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:2 --rg LB:2_AGTACTCATG+GTAGAGTCAG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.AGTACTCATG+GTAGAGTCAG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_2_S2_L002_R1_001.fastq.gz -2 $gbm_data/6931_2_S2_L002_R2_001.fastq.gz -S $gbm_data/alignments/2_rep2.sam

	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:3 --rg LB:3_TACGTGAAGG+GACTGGTTGC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.TACGTGAAGG+GACTGGTTGC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_3_S3_L001_R1_001.fastq.gz -2 $gbm_data/6931_3_S3_L001_R2_001.fastq.gz -S $gbm_data/alignments/3.sam

	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:3 --rg LB:3_TACGTGAAGG+GACTGGTTGC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.TACGTGAAGG+GACTGGTTGC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_3_S3_L002_R1_001.fastq.gz -2 $gbm_data/6931_3_S3_L002_R2_001.fastq.gz -S $gbm_data/alignments/4.sam

	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:4 --rg LB:4_TGTGGTCCGG+GTTCCGCAGG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.TGTGGTCCGG+GTTCCGCAGG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_4_S4_L001_R1_001.fastq.gz -2 $gbm_data/6931_4_S4_L001_R2_001.fastq.gz -S $gbm_data/alignments/5_rep1.sam
	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:4 --rg LB:4_TGTGGTCCGG+GTTCCGCAGG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.TGTGGTCCGG+GTTCCGCAGG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_4_S4_L002_R1_001.fastq.gz -2 $gbm_data/6931_4_S4_L002_R2_001.fastq.gz -S $gbm_data/alignments/5_rep2.sam

	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:5 --rg LB:5_CCGACAGACT+TACCGAACTA --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.CCGACAGACT+TACCGAACTA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_5_S5_L001_R1_001.fastq.gz -2 $gbm_data/6931_5_S5_L001_R2_001.fastq.gz -S $gbm_data/alignments/6_rep1.sam
	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:5 --rg LB:5_CCGACAGACT+TACCGAACTA --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.CCGACAGACT+TACCGAACTA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_5_S5_L002_R1_001.fastq.gz -2 $gbm_data/6931_5_S5_L002_R2_001.fastq.gz -S $gbm_data/alignments/6_rep2.sam

	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:6 --rg LB:6_ATCCAGGTAT+ATCAACAGCC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.ATCCAGGTAT+ATCAACAGCC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_6_S6_L001_R1_001.fastq.gz -2 $gbm_data/6931_6_S6_L001_R2_001.fastq.gz -S $gbm_data/alignments/9.sam

	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:6 --rg LB:6_ATCCAGGTAT+ATCAACAGCC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.ATCCAGGTAT+ATCAACAGCC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_6_S6_L002_R1_001.fastq.gz -2 $gbm_data/6931_6_S6_L002_R2_001.fastq.gz -S $gbm_data/alignments/10.sam



	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:7 --rg LB:7_CGATGCGGTT+GGCTCTTGCG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.CGATGCGGTT+GGCTCTTGCG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_7_S7_L001_R1_001.fastq.gz -2 $gbm_data/6931_7_S7_L001_R2_001.fastq.gz -S $gbm_data/alignments/15.sam

	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:7 --rg LB:7_CGATGCGGTT+GGCTCTTGCG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.CGATGCGGTT+GGCTCTTGCG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_7_S7_L002_R1_001.fastq.gz -2 $gbm_data/6931_7_S7_L002_R2_001.fastq.gz -S $gbm_data/alignments/16.sam


	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:8 --rg LB:8_AATGCGAACA+ATTACTCACC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.AATGCGAACA+ATTACTCACC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_8_S8_L001_R1_001.fastq.gz -2 $gbm_data/6931_8_S8_L001_R2_001.fastq.gz -S $gbm_data/alignments/11.sam

	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:8 --rg LB:8_AATGCGAACA+ATTACTCACC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.AATGCGAACA+ATTACTCACC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_8_S8_L002_R1_001.fastq.gz -2 $gbm_data/6931_8_S8_L002_R2_001.fastq.gz -S $gbm_data/alignments/12.sam



	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:9 --rg LB:9_TTCGGTGTGA+AACAAGGCGT --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.TTCGGTGTGA+AACAAGGCGT -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_9_S9_L001_R1_001.fastq.gz -2 $gbm_data/6931_9_S9_L001_R2_001.fastq.gz -S $gbm_data/alignments/13.sam

	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:9 --rg LB:9_TTCGGTGTGA+AACAAGGCGT --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.TTCGGTGTGA+AACAAGGCGT -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_9_S9_L002_R1_001.fastq.gz -2 $gbm_data/6931_9_S9_L002_R2_001.fastq.gz -S $gbm_data/alignments/14.sam


	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:10 --rg LB:10_CATTAACTGA+AGAATCTTCG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.CATTAACTGA+AGAATCTTCG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_10_S10_L001_R1_001.fastq.gz -2 $gbm_data/6931_10_S10_L001_R2_001.fastq.gz -S $gbm_data/alignments/17.sam

	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:10 --rg LB:10_CATTAACTGA+AGAATCTTCG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.CATTAACTGA+AGAATCTTCG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_10_S10_L002_R1_001.fastq.gz -2 $gbm_data/6931_10_S10_L002_R2_001.fastq.gz -S $gbm_data/alignments/18.sam


	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:11 --rg LB:11_CCGGTTCCTA+CGGCAATGGA --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.CCGGTTCCTA+CGGCAATGGA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_11_S11_L001_R1_001.fastq.gz -2 $gbm_data/6931_11_S11_L001_R2_001.fastq.gz -S $gbm_data/alignments/19.sam

	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:11 --rg LB:11_CCGGTTCCTA+CGGCAATGGA --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.CCGGTTCCTA+CGGCAATGGA -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_11_S11_L002_R1_001.fastq.gz -2 $gbm_data/6931_11_S11_L002_R2_001.fastq.gz -S $gbm_data/alignments/20.sam


	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:12 --rg LB:12_ATACATCACA+TAACAGTGTT --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.ATACATCACA+TAACAGTGTT -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_12_S12_L001_R1_001.fastq.gz -2 $gbm_data/6931_12_S12_L001_R2_001.fastq.gz -S $gbm_data/alignments/21_rep1.sam
	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:12 --rg LB:12_ATACATCACA+TAACAGTGTT --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.ATACATCACA+TAACAGTGTT -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_12_S12_L002_R1_001.fastq.gz -2 $gbm_data/6931_12_S12_L002_R2_001.fastq.gz -S $gbm_data/alignments/21_rep2.sam


	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:13 --rg LB:13_TAGAGAATAC+AGAGTGCGGC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.TAGAGAATAC+AGAGTGCGGC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_13_S13_L001_R1_001.fastq.gz -2 $gbm_data/6931_13_S13_L001_R2_001.fastq.gz -S $gbm_data/alignments/22_rep1.sam
	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:13 --rg LB:13_TAGAGAATAC+AGAGTGCGGC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.TAGAGAATAC+AGAGTGCGGC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_13_S13_L002_R1_001.fastq.gz -2 $gbm_data/6931_13_S13_L002_R2_001.fastq.gz -S $gbm_data/alignments/22_rep2.sam


	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:14 --rg LB:14_GTCTCGCCAC+CAGCCGATTG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.GTCTCGCCAC+CAGCCGATTG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_14_S14_L001_R1_001.fastq.gz -2 $gbm_data/6931_14_S14_L001_R2_001.fastq.gz -S $gbm_data/alignments/23.sam

	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:14 --rg LB:14_GTCTCGCCAC+CAGCCGATTG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.GTCTCGCCAC+CAGCCGATTG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_14_S14_L002_R1_001.fastq.gz -2 $gbm_data/6931_14_S14_L002_R2_001.fastq.gz -S $gbm_data/alignments/24.sam


	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:15 --rg LB:15_CGCTGTCTCA+ATGTCGTGGT --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.CGCTGTCTCA+ATGTCGTGGT -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_15_S15_L001_R1_001.fastq.gz -2 $gbm_data/6931_15_S15_L001_R2_001.fastq.gz -S $gbm_data/alignments/25.sam

	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:15 --rg LB:15_CGCTGTCTCA+ATGTCGTGGT --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.CGCTGTCTCA+ATGTCGTGGT -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_15_S15_L002_R1_001.fastq.gz -2 $gbm_data/6931_15_S15_L002_R2_001.fastq.gz -S $gbm_data/alignments/26.sam


	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:16 --rg LB:16_GCTAATAGGA+AGAGCACTAG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.GCTAATAGGA+AGAGCACTAG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_16_S16_L001_R1_001.fastq.gz -2 $gbm_data/6931_16_S16_L001_R2_001.fastq.gz -S $gbm_data/alignments/27.sam

	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:16 --rg LB:16_GCTAATAGGA+AGAGCACTAG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.GCTAATAGGA+AGAGCACTAG -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_16_S16_L002_R1_001.fastq.gz -2 $gbm_data/6931_16_S16_L002_R2_001.fastq.gz -S $gbm_data/alignments/28.sam


	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:20 --rg LB:20_GCGTGCTGTG+CGGTGACACC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.GCGTGCTGTG+CGGTGACACC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_20_S20_L001_R1_001.fastq.gz -2 $gbm_data/6931_20_S20_L001_R2_001.fastq.gz -S $gbm_data/alignments/7_rep1.sam
	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:20 --rg LB:20_GCGTGCTGTG+CGGTGACACC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.GCGTGCTGTG+CGGTGACACC -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_20_S20_L002_R1_001.fastq.gz -2 $gbm_data/6931_20_S20_L002_R2_001.fastq.gz -S $gbm_data/alignments/7_rep2.sam


	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:21 --rg LB:21_CGAGAGGCGT+GAGACATAAT --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.CGAGAGGCGT+GAGACATAAT -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_21_S21_L001_R1_001.fastq.gz -2 $gbm_data/6931_21_S21_L001_R2_001.fastq.gz -S $gbm_data/alignments/8_rep1.sam
	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:21 --rg LB:21_CGAGAGGCGT+GAGACATAAT --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.CGAGAGGCGT+GAGACATAAT -x $gbm/RNA_REF_FA/hg38/genome --dta --rna-strandness FR -1 $gbm_data/6931_21_S21_L002_R1_001.fastq.gz -2 $gbm_data/6931_21_S21_L002_R2_001.fastq.gz -S $gbm_data/alignments/8_rep2.sam


the following samples are E0771 mouse breast cancer mets to the brain, so need to align with mouse reference genome


	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:17 --rg LB:17_CCTAACACAG+GGTGGAATAC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.CCTAACACAG+GGTGGAATAC -x $gbm/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $gbm_data/6931_17_S17_L001_R1_001.fastq.gz -2 $gbm_data/6931_17_S17_L001_R2_001.fastq.gz -S $gbm_data/alignments/31_rep1.sam
	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:17 --rg LB:17_CCTAACACAG+GGTGGAATAC --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.CCTAACACAG+GGTGGAATAC -x $gbm/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $gbm_data/6931_17_S17_L002_R1_001.fastq.gz -2 $gbm_data/6931_17_S17_L002_R2_001.fastq.gz -S $gbm_data/alignments/31_rep2.sam


	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:18 --rg LB:18_TGCCGGTCAG+GAGGCTCCTA --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.TGCCGGTCAG+GAGGCTCCTA -x $gbm/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $gbm_data/6931_18_S18_L001_R1_001.fastq.gz -2 $gbm_data/6931_18_S18_L001_R2_001.fastq.gz -S $gbm_data/alignments/32_rep1.sam
	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:18 --rg LB:18_TGCCGGTCAG+GAGGCTCCTA --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.TGCCGGTCAG+GAGGCTCCTA -x $gbm/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $gbm_data/6931_18_S18_L002_R1_001.fastq.gz -2 $gbm_data/6931_18_S18_L002_R2_001.fastq.gz -S $gbm_data/alignments/32_rep2.sam


	hisat2 -p 8 --rg-id=HCTMMDRXY.1 --rg SM:19 --rg LB:19_TTAACCTTCG+TAATGGCAAG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.1.TTAACCTTCG+TAATGGCAAG -x $gbm/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $gbm_data/6931_19_S19_L001_R1_001.fastq.gz -2 $gbm_data/6931_19_S19_L001_R2_001.fastq.gz -S $gbm_data/alignments/29.sam

	hisat2 -p 8 --rg-id=HCTMMDRXY.2 --rg SM:19 --rg LB:19_TTAACCTTCG+TAATGGCAAG --rg PL:ILLUMINA --rg PU:HCTMMDRXY.2.TTAACCTTCG+TAATGGCAAG -x $gbm/RNA_REF_FA/mm10/genome --dta --rna-strandness FR -1 $gbm_data/6931_19_S19_L002_R1_001.fastq.gz -2 $gbm_data/6931_19_S19_L002_R2_001.fastq.gz -S $gbm_data/alignments/30.sam




   


## Convert SAM files into BAM files

	cd $gbm_data/alignments

	samtools sort -@ 8 -o 1_rep1.bam 1_rep1.sam
	samtools sort -@ 8 -o 1_rep2.bam 1_rep2.sam

	samtools sort -@ 8 -o 2_rep1.bam 2_rep1.sam
	samtools sort -@ 8 -o 2_rep2.bam 2_rep2.sam

	samtools sort -@ 8 -o 3.bam 3.sam

	samtools sort -@ 8 -o 4.bam 4.sam

	samtools sort -@ 8 -o 5_rep1.bam 5_rep1.sam
	samtools sort -@ 8 -o 5_rep2.bam 5_rep2.sam

	samtools sort -@ 8 -o 6_rep1.bam 6_rep1.sam
	samtools sort -@ 8 -o 6_rep2.bam 6_rep2.sam

	samtools sort -@ 8 -o 7_rep1.bam 7_rep1.sam
	samtools sort -@ 8 -o 7_rep2.bam 7_rep2.sam

	samtools sort -@ 8 -o 8_rep1.bam 8_rep1.sam
	samtools sort -@ 8 -o 8_rep2.bam 8_rep2.sam





	samtools sort -@ 8 -o 9.bam 9.sam

	samtools sort -@ 8 -o 10.bam 10.sam

	samtools sort -@ 8 -o 11.bam 11.sam

	samtools sort -@ 8 -o 12.bam 12.sam

	samtools sort -@ 8 -o 13.bam 13.sam

	samtools sort -@ 8 -o 14.bam 14.sam

	samtools sort -@ 8 -o 15.bam 15.sam

	samtools sort -@ 8 -o 16.bam 16.sam




	samtools sort -@ 8 -o 17.bam 17.sam

	samtools sort -@ 8 -o 18.bam 18.sam

	samtools sort -@ 8 -o 19.bam 19.sam

	samtools sort -@ 8 -o 20.bam 20.sam

	samtools sort -@ 8 -o 21_rep1.bam 21_rep1.sam
	samtools sort -@ 8 -o 21_rep2.bam 21_rep2.sam

	samtools sort -@ 8 -o 22_rep1.bam 22_rep1.sam
	samtools sort -@ 8 -o 22_rep2.bam 22_rep2.sam





	samtools sort -@ 8 -o 23.bam 23.sam

	samtools sort -@ 8 -o 24.bam 24.sam

	samtools sort -@ 8 -o 25.bam 25.sam

	samtools sort -@ 8 -o 26.bam 26.sam

	samtools sort -@ 8 -o 27.bam 27.sam

	samtools sort -@ 8 -o 28.bam 28.sam




	samtools sort -@ 8 -o 29.bam 29.sam

	samtools sort -@ 8 -o 30.bam 30.sam

	samtools sort -@ 8 -o 31_rep1.bam 31_rep1.sam
	samtools sort -@ 8 -o 31_rep2.bam 31_rep2.sam

	samtools sort -@ 8 -o 32_rep1.bam 32_rep1.sam
	samtools sort -@ 8 -o 32_rep2.bam 32_rep2.sam



## Merge the bam files


	java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=1.bam INPUT=1_rep1.bam INPUT=1_rep2.bam
	java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=2.bam INPUT=2_rep1.bam INPUT=2_rep2.bam

	java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=5.bam INPUT=5_rep1.bam INPUT=5_rep2.bam
	java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=6.bam INPUT=6_rep1.bam INPUT=6_rep2.bam
	java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=7.bam INPUT=7_rep1.bam INPUT=7_rep2.bam
	java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=8.bam INPUT=8_rep1.bam INPUT=8_rep2.bam

	java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=21.bam INPUT=21_rep1.bam INPUT=21_rep2.bam
	java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=22.bam INPUT=22_rep1.bam INPUT=22_rep2.bam

	java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=31.bam INPUT=31_rep1.bam INPUT=31_rep2.bam
	java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=32.bam INPUT=32_rep1.bam INPUT=32_rep2.bam







## Index all bam files for IGV visualization
Make sure that the only files in the directory are the sam and bam/bai files, then:

can do:
samtools index 1.bam

or:

    find *.bam -exec echo samtools index {} \; | sh



## use samtoools flagstat to geta basic sumary of an alignment 
for example percent of unmapped reads

	mkdir flagstat
	samtools flagstat 1.bam > flagstat/1.bam.flagstat
	samtools flagstat 2.bam > flagstat/2.bam.flagstat


view the resulting flagstat: 
	
	cat flagstat/1.bam.flagstat
	


## Use Stringtie to generate expression estimates (FPKM and TPM) from the SAM/BAM files generated by HISAT2 and samtools
this takes the bam file, add expression estimates FPKM and TPM using the reference annotation GTF file as a guide, and gives an annotated gtf file for the sample.
    
  
	mkdir -p $gbm/expression/stringtie/ref_only/
	cd $gbm/expression/stringtie/ref_only/


	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 1/transcripts.gtf -A 1/gene_abundances.tsv $gbm/alignments/1.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 2/transcripts.gtf -A 2/gene_abundances.tsv $gbm/alignments/2.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 3/transcripts.gtf -A 3/gene_abundances.tsv $gbm/alignments/3.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 4/transcripts.gtf -A 4/gene_abundances.tsv $gbm/alignments/4.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 5/transcripts.gtf -A 5/gene_abundances.tsv $gbm/alignments/5.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 6/transcripts.gtf -A 6/gene_abundances.tsv $gbm/alignments/6.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 7/transcripts.gtf -A 7/gene_abundances.tsv $gbm/alignments/7.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 8/transcripts.gtf -A 8/gene_abundances.tsv $gbm/alignments/8.bam

	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 9/transcripts.gtf -A 9/gene_abundances.tsv $gbm/alignments/9.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 10/transcripts.gtf -A 10/gene_abundances.tsv $gbm/alignments/10.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 11/transcripts.gtf -A 11/gene_abundances.tsv $gbm/alignments/11.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 12/transcripts.gtf -A 12/gene_abundances.tsv $gbm/alignments/12.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 13/transcripts.gtf -A 13/gene_abundances.tsv $gbm/alignments/13.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 14/transcripts.gtf -A 14/gene_abundances.tsv $gbm/alignments/14.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 15/transcripts.gtf -A 15/gene_abundances.tsv $gbm/alignments/15.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 16/transcripts.gtf -A 16/gene_abundances.tsv $gbm/alignments/16.bam


	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 17/transcripts.gtf -A 17/gene_abundances.tsv $gbm/alignments/17.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 18/transcripts.gtf -A 18/gene_abundances.tsv $gbm/alignments/18.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 19/transcripts.gtf -A 19/gene_abundances.tsv $gbm/alignments/19.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 20/transcripts.gtf -A 20/gene_abundances.tsv $gbm/alignments/20.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 21/transcripts.gtf -A 21/gene_abundances.tsv $gbm/alignments/21.bam

	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 22/transcripts.gtf -A 22/gene_abundances.tsv $gbm/alignments/22.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 23/transcripts.gtf -A 23/gene_abundances.tsv $gbm/alignments/23.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 24/transcripts.gtf -A 24/gene_abundances.tsv $gbm/alignments/24.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 25/transcripts.gtf -A 25/gene_abundances.tsv $gbm/alignments/25.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 26/transcripts.gtf -A 26/gene_abundances.tsv $gbm/alignments/26.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 27/transcripts.gtf -A 27/gene_abundances.tsv $gbm/alignments/27.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf -e -B -o 28/transcripts.gtf -A 28/gene_abundances.tsv $gbm/alignments/28.bam



	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/mm10.ncbiRefSeq.gtf -e -B -o 29/transcripts.gtf -A 29/gene_abundances.tsv $gbm/alignments/29.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/mm10.ncbiRefSeq.gtf -e -B -o 30/transcripts.gtf -A 30/gene_abundances.tsv $gbm/alignments/30.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/mm10.ncbiRefSeq.gtf -e -B -o 31/transcripts.gtf -A 31/gene_abundances.tsv $gbm/alignments/31.bam
	stringtie --fr -p 8 -G $gbm/RNA_REF_GTF/mm10.ncbiRefSeq.gtf -e -B -o 32/transcripts.gtf -A 32/gene_abundances.tsv $gbm/alignments/32.bam





view raw output from stringtie

	less -S 1/transcripts.gtf
	
View transcript records only (column -t makes formatting easier to view). scroll right to see the ### 3 major expression matrices "cov, fpkm, and tpm". ### Cov, or coverage, is average per base covereage over the genomic segment. basically, the coverage of each transcript. it's un-normalized unlike fpkm or tpm.

	grep -v "^#" 1/transcripts.gtf | grep -w "transcript" | column -t | less -S

Limit the view only to transcript records and their expression estimates (FPKM and TPM values). the if statement {if ($3=="transcript") print} means only print the line if the third column is a transcript. then only return field 1, 4, 9 which are the fpkm nad tpm values.

	awk '{if ($3=="transcript") print}' 1/transcripts.gtf | cut -f 1,4,9 | less -S

Could also view gene and trnscript level expression values in the two files generated by stringtie, a transcript.gtf file (transcript level abundance) and a gene_abundances.tsv (a gene level abundance):

	column -t 1/t_data.ctab | less -S
	less -S -x20 1/gene_abundances.tsv


# htseq for raw counts DE
## Parallel to Stringtie (produces FPKM/TPM), can also run htseq-count on hisat2 alignments to produce raw counts (instead of FPKM/TPM values) for differential expression analysis


htseq-count basic usage:

	htseq-count [options] <sam_file> <gff_file>

Extra options specified below:

	’–format’ specify the input file format one of BAM or SAM. Since we have BAM format files, select ‘bam’ for this option.
	’–order’ provide the expected sort order of the input file. Previously we generated position sorted BAM files so use ‘pos’.
	’–mode’ determines how to deal with reads that overlap more than one feature. We believe the ‘intersection-strict’ mode is best.
	’–stranded’ specifies whether data is stranded or not. The TruSeq strand-specific RNA libraries suggest the ‘reverse’ option for this parameter.
	’–minaqual’ will skip all reads with alignment quality lower than the given minimum value
	’–type’ specifies the feature type (3rd column in GFF file) to be used. (default, suitable for RNA-Seq and Ensembl GTF files: exon)
	’–idattr’ The feature ID used to identify the counts in the output table. The default, suitable for RNA-SEq and Ensembl GTF files, is gene_id.


	mkdir -p $gbm_data/expression/htseq_counts
	cd $gbm_data/expression/htseq_counts


	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/1.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 1.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/2.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 2.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/3.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 3.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/4.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 4.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/5.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 5.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/6.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 6.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/7.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 7.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/8.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 8.tsv

	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/9.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 9.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/10.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 10.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/11.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 11.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/12.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 12.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/13.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 13.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/14.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 14.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/15.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 15.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/16.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 16.tsv


	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/17.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 17.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/18.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 18.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/19.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 19.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/20.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 20.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/21.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 21.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/22.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 22.tsv


	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/23.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 23.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/24.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 24.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/25.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 25.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/26.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 26.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/27.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 27.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/28.bam $gbm/RNA_REF_GTF/hg38.ncbiRefSeq.gtf > 28.tsv


	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/29.bam $gbm/RNA_REF_GTF/mm10.ncbiRefSeq.gtf > 29.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/30.bam $gbm/RNA_REF_GTF/mm10.ncbiRefSeq.gtf > 30.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/31.bam $gbm/RNA_REF_GTF/mm10.ncbiRefSeq.gtf > 31.tsv
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $gbm_data/alignments/32.bam $gbm/RNA_REF_GTF/mm10.ncbiRefSeq.gtf > 32.tsv


Join the results for each replicate together and merge results files into a single matrix for use in edgeR. Since edgeR only does pairwise comparisons, need to pair up the compared samples, then creat a simple text file with just the header that will be used for the table:


slice vs invitro: 

	join 1.tsv 2.tsv | join - 7.tsv | join - 8.tsv > GBM011_slice_vs_invitro_gene_read_counts_table_all.tsv
	echo "GeneID 1 2 7 8" > 011_slice_vs_invitro_header.txt

Clean up a bit more, add a header, reformat the result as a tab delimited file. note: grep -v "__" is being used to filter out the summary lines at the end of the files that ht-seq count gives to summarize reads that had no feature, were ambiguous, did not align at all, did not align due to poor alignment quality, or the alignment was not unique.

(note: awk -v OFS="\t" '$1=$1' is using awk to replace the single space characters that were in the concatenated version of our header.txt and gene_read_counts_table_all.tsv with a tab character. -v is used to reset the variable OFS, which stands for Output Field Separator. By default, this is a single space. By specifying OFS="\t", we are telling awk to replace the single space with a tab. The '$1=$1' tells awk to reevaluate the input using the new output variable)


	cat GBM011_slice_vs_invitro_header.txt GBM011_slice_vs_invitro_gene_read_counts_table_all.tsv | grep -v "__" | awk -v OFS="\t" '$1=$1' > 011_slice_vs_invitro_gene_read_counts_table_all_final.tsv


	head GBM011_slice_vs_invitro_gene_read_counts_table_all_final.tsv | column -t

	rm -f 011_slice_vs_invitro_header.txt 011_slice_vs_invitro_header.txt




tissue vs invitro:

organoid vs invitro:




invitro vs slice:

tissue vs slice:

organoid vs slice:




slice vs tissue:

organoid vs tissue:

invitro vs tissue:





copy all of above for 024 sample:



(DE of htseq results using EdgeR continues after Stringtie+Ballgown)


-------------------------



# Use Ballgown in R for differential expression (DE) analysis (then PCA) using output from Stringtie
## Perform A vs. B comparison, using all replicates, for known (reference only mode) transcripts

(raw counts using htseq output will come later)

	mkdir -p $gbm/de/ballgown/ref_only
	cd $gbm/de/ballgown/ref_only/


Use printf to create/print a table with ids, type (each type of sample is a type), and path to the file, as the header. Then n returns a new line. 
## (note: id needs to match the file folder names created by stringtie)

Bascially, need a table that needs to look like this to feed into R:

ids type path-to-file-011_invitro_1 011 $gbm/expression/stringtie/1 011_invitro_2 011 $gbm/expression/stringtie/2 ... ...

goal is to generate a header file to load into R, for ALL samples for principal component analysis (the simplest form of multidimentional scaling), and also a file for pairwise comparisons. since we have a ton of comparisisons, might just not do this for now and only do the PCA. 

file for all 011 samples for PCA: (this is how the script should look like (without the enters inbetween each line):

printf "\"ids\",\"type\",\"path
\"\n\"1\",\"011_slice\",\"$gbm/expression/stringtie/ref_only/1
\"\n\"2\",\"011_slice\",\"$gbm/expression/stringtie/ref_only/2
\"\n\"3\",\"011_organoid\",\"$gbm/expression/stringtie/ref_only/3
\"\n\"4\",\"011_organoid\",\"$gbm/expression/stringtie/ref_only/4
\"\n\"5\",\"011_tissue\",\"$gbm/expression/stringtie/ref_only/5
\"\n\"6\",\"011_tissue\",\"$gbm/expression/stringtie/ref_only/6
\"\n\"7\",\"011_invitro\",\"$gbm/expression/stringtie/ref_only/7
\"\n\"8\",\"011_invitro\",\"$gbm/expression/stringtie/ref_only/8
\"\n" > GBM011_all.csv

	cd $gbm/de/ballgown/ref_only/

	printf "\"ids\",\"type\",\"path\"\n\"1\",\"011_slice\",\"$gbm/expression/stringtie/ref_only/1\"\n\"2\",\"011_slice\",\"$gbm/expression/stringtie/ref_only/2\"\n\"3\",\"011_organoid\",\"$gbm/expression/stringtie/ref_only/3\"\n\"4\",\"011_organoid\",\"$gbm/expression/stringtie/ref_only/4\"\n\"5\",\"011_tissue\",\"$gbm/expression/stringtie/ref_only/5\"\n\"6\",\"011_tissue\",\"$gbm/expression/stringtie/ref_only/6\"\n\"7\",\"011_invitro\",\"$gbm/expression/stringtie/ref_only/7\"\n\"8\",\"011_invitro\",\"$gbm/expression/stringtie/ref_only/8\"\n" > GBM011_all.csv

R script:


	R --no-restore
	library(ballgown)
	library(genefilter)
	library(dplyr)
	library(devtools)
	library(ggplot2)
	library(gplots)
	library(GenomicRanges)

	pheno_data = read.csv("GBM011_all.csv")  


	bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)
	bg


	bg_table = texpr(bg, 'all')


	bg_gene_names = unique(bg_table[, 9:10])
	head(bg_gene_names)


	save(bg, file='bg.rda')

	bg
output: ballgown instance with 190734 transcripts and 8 samples


	pdf(file="GBM011_R_output.pdf")

	working_dir = "~/workspace/gbm/de/ballgown/ref_only"
	setwd(working_dir)
	dir()


Import expression and differential expression results from the HISAT2/StringTie/Ballgown pipeline
	
	load('bg.rda')


Load gene names for lookup later in the tutorial
	
	bg_table = texpr(bg, 'all')
	bg_gene_names = unique(bg_table[, 9:10])


Pull the gene_expression data frame from the ballgown object

	gene_expression = as.data.frame(gexpr(bg))
	head(gene_expression)

View the column names

	colnames(gene_expression)

View the row names
	
	row.names(gene_expression)

Determine the dimensions of the dataframe.  'dim()' will return the number of rows and columns

	dim(gene_expression)


Just for fun, check BRD4 expression across all 8 samples:

	i = row.names(gene_expression) == "BRD4"
	gene_expression[i,]






Load the transcript to gene index from the ballgown object. Each row of data represents a transcript. Many of these transcripts represent the same gene. Determine the numbers of transcripts and unique genes  


	transcript_gene_table = indexes(bg)$t2g
	head(transcript_gene_table)
	
	length(row.names(transcript_gene_table)) #Transcript count
	length(unique(transcript_gene_table[,"g_id"])) #Unique Gene count

> length(row.names(transcript_gene_table)) #Transcript count
[1] 190734
> length(unique(transcript_gene_table[,"g_id"])) #Unique Gene count
[1] 54651




Plot the number of transcripts per gene. Many genes will have only 1 transcript, some genes will have several transcripts. Use the 'table()' command to count the number of times each gene symbol occurs (i.e. the # of transcripts that have each gene symbol). Then use the 'hist' command to create a histogram of these counts

	counts=table(transcript_gene_table[,"g_id"])
	c_one = length(which(counts == 1))
	c_more_than_one = length(which(counts > 1))
	c_max = max(counts)
	hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
	legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
	legend("topright", legend_text, lty=NULL)


Plot the distribution of transcript sizes as a histogram. lengths will be those of known transcripts. Good QC step: we had a low coverage library, or other problems, we might get short 'transcripts' that are actually only pieces of real transcripts.

	full_table <- texpr(bg , 'all')
	hist(full_table$length, breaks=500, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")

View the summary FPKM values (minimum and maximum FPKM values) for any particular library

	min(gene_expression[,"FPKM.1"])
	max(gene_expression[,"FPKM.2"])


Set the minimum non-zero FPKM values by one of two ways:

coverting 0's to NA, and calculating the minimum or all non NA values
two ways: 
zz = fpkm_matrix[,data_columns]
zz[zz==0] = NA
min_nonzero = min(zz, na.rm=TRUE)
min_nonzero


Alternatively just set min value to 1
	
	min_nonzero=1

Set the columns for finding FPKM and create shorter names for figures

	data_columns=c(1:8)
	short_names=c("slice_1","slice2","organoid_1","organoid_2","tissue_1","tissue_2","invitro_1","invitro_2")


Plot range of values and general distribution of FPKM values for all 8 libraries

Create boxplots using different colors by setting storing the colors of the columns in a variable called data_colors. then display on a log2 scale and add the minimum non-zero value to avoid log2(0). Note that the bold horizontal line on each boxplot is the median.


	colors()
	data_colors=c("tomato1","tomato2","royalblue1","royalblue2","seagreen1","seagreen2","grey1","grey2")

	boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all 8 sample libraries")



## plot a pair of replicates to assess reproducibility of technical replicates. 
Tranform the data by converting to log2 scale after adding an arbitrary small value to avoid log2(0). Also add a straight line of slope 1, and intercept 0. Also calculate the correlation coefficient and display in a legend.

	x = gene_expression[,"FPKM.1"]
	y = gene_expression[,"FPKM.2"]
	plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (011_slice, Replicate 1)", ylab="FPKM (011_slice, Replicate 2)", main="Comparison of expression values for replicates of slice samples")
	abline(a=0,b=1)
	rs=cor(x,y)^2
	legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")


check organoid samples

	x = gene_expression[,"FPKM.3"]
	y = gene_expression[,"FPKM.4"]
	plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (011_organoid, Replicate 1)", ylab="FPKM (011_organoid, Replicate 2)", main="Comparison of expression values for replicates of organoid samples")
	abline(a=0,b=1)
	rs=cor(x,y)^2
	legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")

check tissue samples

	x = gene_expression[,"FPKM.5"]
	y = gene_expression[,"FPKM.6"]
	plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (011_tissue, Replicate 1)", ylab="FPKM (011_tissue, Replicate 2)", main="Comparison of expression values for replicates of tissue samples")
	abline(a=0,b=1)
	rs=cor(x,y)^2
	legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")


check in vitro samples

	x = gene_expression[,"FPKM.7"]
	y = gene_expression[,"FPKM.8"]
	plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (011_invitro, Replicate 1)", ylab="FPKM (011_invitro, Replicate 2)", main="Comparison of expression values for replicates of in vitro samples")
	abline(a=0,b=1)
	rs=cor(x,y)^2
	legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")



## Compare the correlation distance between all replicates

Calculate the FPKM sum for all 8 libraries

	gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)

Filter out genes with a grand sum FPKM of less than 10

	i = which(gene_expression[,"sum"] > 10)


Calculate the correlation between all pairs of data

	r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")
	r
	
	
## Plot MDS.
Convert correlation to distance, and use 'multi-dimensional scaling' to plot the relative differences between libraries, by calculating 2-dimensional coordinates to plot points for each library using eigenvectors (eig=TRUE). d, k=2 means 2 dimensions
	
	d=1-r
	mds=cmdscale(d, k=2, eig=TRUE)
	par(mfrow=c(1,1))
	plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes)", xlim=c(-0.4,0.4), ylim=c(-0.4,0.4))
	points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
	text(mds$points[,1], mds$points[,2], short_names, col=data_colors)

Calculate the differential expression results including significance

	results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
	results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))



close out the PDF
dev.off()



R script for other samples are located in other files in the same danluo7/gbm project folder




# Could also do pairwise comparisons for 011 (but has to be two at a time):

against invitro:
printf "\"ids\",\"type\",\"path
\"\n\"1\",\"011_slice\",\"$gbm/expression/stringtie/1
\"\n\"2\",\"011_slice\",\"$gbm/expression/stringtie/2
\"\n\"7\",\"011_invitro\",\"$gbm/expression/stringtie/7
\"\n\"8\",\"011_invitro\",\"$gbm/expression/stringtie/8
\"\n" > 011_slice_vs_invitro.csv

	printf "\"ids\",\"type\",\"path\"\n\"1\",\"011_slice\",\"$gbm/expression/stringtie/1\"\n\"2\",\"011_slice\",\"$gbm/expression/stringtie/2\"\n\"7\",\"011_invitro\",\"$gbm/expression/stringtie/7\"\n\"8\",\"011_invitro\",\"$gbm/expression/stringtie/8\"\n" > 011_slice_vs_invitro.csv

	
R script (for slice vs in vitro comparison), rest in gbm folder.

	R --no-restore
	library(ballgown)
	library(genefilter)
	library(dplyr)
	library(devtools)

	pheno_data = read.csv("011_slice_vs_invitro.csv")     ###change me


	bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)
	bg


	bg_table = texpr(bg, 'all')


	bg_gene_names = unique(bg_table[, 9:10])
	head(bg_gene_names)


	save(bg, file='bg.rda')
	
The stattest function below will require pairwise comparisons and can't know multiple comparisons at once

	results_transcripts = stattest(bg, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")


	results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")

	head(results_genes)


	results_genes = merge(results_genes, bg_gene_names, by.x=c("id"), by.y=c("gene_id"))

	write.table(results_transcripts, "011_slice_vs_invitro_transcript_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)    ###change me
	write.table(results_genes, "011_slice_vs_invitro_gene_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)                 ###change me


	bg_filt = subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)


	bg_filt_table = texpr(bg_filt , 'all')
	bg_filt_gene_names = unique(bg_filt_table[, 9:10])


	results_transcripts = stattest(bg_filt, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
	results_genes = stattest(bg_filt, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
	results_genes = merge(results_genes, bg_filt_gene_names, by.x=c("id"), by.y=c("gene_id"))


	write.table(results_transcripts, "011_slice_vs_invitro_transcript_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)   ###change me
	write.table(results_genes, "011_slice_vs_invitro_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)    ###change me
	
	sig_transcripts = subset(results_transcripts, results_transcripts$pval<0.05)
	sig_genes = subset(results_genes, results_genes$pval<0.05)

	head(sig_genes)

	nrow(sig_genes)

	write.table(sig_transcripts, "011_slice_vs_invitro_transcript_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)      ###change me
	write.table(sig_genes, "011_slice_vs_invitro_gene_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)         ###change me

	quit()
	n

	grep -v feature 011_slice_vs_invitro_gene_results_filtered.tsv | wc -l

	grep -v feature 011_slice_vs_invitro_gene_results_sig.tsv | sort -rnk 3 | head -n 20 | column -t 	#Higher abundance in invitro
	grep -v feature 011_slice_vs_invitro_gene_results_sig.tsv | sort -nk 3 | head -n 20 | column -t 	#Higher abundance in tissue



	grep -v feature 011_slice_vs_invitro_gene_results_sig.tsv | cut -f 6 | sed 's/\"//g' > 011_slice_vs_invitro_DE_genes.txt

	head 011_slice_vs_invitro_DE_genes.txt

11 more scripts in gbm folder.









# Parallel to Ballgown, also need to use edgeR for DE analysis

	cd $gbm
	mkdir -p de/htseq_counts
	cd de/htseq_counts

## Launch R, set working directory, read in the count matrix tsv file created by htseq, and check dimensions of the file.

	R --no-restore
	working_dir = "~/workspace/gbm/de/htseq_counts"
	setwd(working_dir)


Doing 011 samples for slice vs in vitro pairwise comparisons first.

	rawdata=read.table("~/workspace/gbm/expression/htseq_counts/GBM011_slice_vs_invitro_gene_read_counts_table_all_final.tsv", header=TRUE, stringsAsFactors=FALSE, row.names=1)

	dim(rawdata)

output: [1] 54651     4


Require at least 1/4 of samples to have expressed count >= 10

	sample_cutoff <- (1/4)
	count_cutoff <- 10

Define a function to calculate the fraction of values expressed above the count cutoff

	getFE <- function(data,count_cutoff){
	 FE <- (sum(data >= count_cutoff)/length(data))
	 return(FE)
	}

Apply the function to all genes, and filter out genes not meeting the sample cutoff

	fraction_expressed <- apply(rawdata, 1, getFE, count_cutoff)
	keep <- which(fraction_expressed >= sample_cutoff)
	rawdata <- rawdata[keep,]

Check dimensions again to see effect of filtering

	dim(rawdata)

output: 

load edgeR
[1] 18302     4
went from 54k to 18k. 

	library('edgeR')

make class labels

	class <- c( rep("011_invitro",2), rep("011_slice",2) )


Get common gene names (placeholder code for when I need to use ENSEMBL reference genomes instead of NCBI) in that case, will just need to do Gene=rownames(rawdata) Symbol=mapping[Gene,1] gene_annotations=cbind(Gene,Symbol)

	Gene=rownames(rawdata)
	Symbol=rownames(rawdata)
	gene_annotations=cbind(Gene,1)

Make DGEList object

	y <- DGEList(counts=rawdata, genes=gene_annotations, group=class)
	nrow(y)

TMM Normalization

	y <- calcNormFactors(y)

Estimate dispersion

	y <- estimateCommonDisp(y, verbose=TRUE)
	y <- estimateTagwiseDisp(y)

Differential expression test

	et <- exactTest(y)

Extract raw counts to add back onto DE results

	counts <- getCounts(y)

Print top genes

	topTags(et)

Print number of up/down significant genes at FDR = 0.05 significance level

	summary(de <- decideTestsDGE(et, adjust.method="BH", p=.05))

Get output with BH-adjusted FDR values - all genes, any p-value, unsorted

	out <- topTags(et, n = "Inf", adjust.method="BH", sort.by="none", p.value=1)$table

Add raw counts back onto results for convenience (make sure sort and total number of elements allows proper join)

	out2 <- cbind(out, counts)

Limit to significantly DE genes

	out3 <- out2[as.logical(de),]

Order by q-value

	o <- order(et$table$PValue[as.logical(de)],decreasing=FALSE)

	out4 <- out3[o,]

Save table

	write.table(out4, file="011_DE_genes.txt", quote=FALSE, row.names=FALSE, sep="\t")

q() then

	cd $gbm/de/htseq_counts/

	cut -f 1 $gbm/de/htseq_counts/011_DE_genes.txt | sort | uniq | grep -v Gene_Name > 011_DE_gene_symbols.txt



(not sure what happened to the video of first half of this tutorial: https://rnabio.org/module-03-expression/0003/04/01/DE_Visualization/)






starting from "supplementary R analysis:

# creat multidimention plots to visualize differences between samples and replicates within samples

	cd $gbm/de/ballgown/ref_only/

R scripts in same folder.


