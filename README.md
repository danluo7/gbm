# gbm

## check quality of reads

    fastqc *.fastq.gz 
    
    python3 -m multiqc . 
    
    mkdir fastqc
    mv *fastqc* fastqc
    
Reference genome:    
  
Sources for obtaining gene annotation files formatted for HISAT2: HISAT2 Precomputed Genome Index (used by Andrew in Hingtgen's lab) available from their FTP site ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/

    wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/hg38.tar.gz
    tar -xzvf hg38.tar.gz
  
These 8 files together constitute the index: they are all that is needed to align reads to that reference. 

Download the reference genome from UCSC: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/. File used by Andrew is: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz

    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz 
    gzip -d hg38.ncbiRefSeq.gtf.gz


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


example: 
hisat2 -p 8 --rg-id=UHR_Rep1 --rg SM:UHR --rg LB:UHR_Rep1_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-ACTGAC.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz -S ./UHR_Rep1.sam






