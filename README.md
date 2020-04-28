# ExACvcf:

ExACvcf is a prototype variant annotation tool. It annotates VCF files through Broad Institute ExAC Project and output a new VCF file annotating each variant in the file with additional information from ExAC. Each variant is annotated with the following pieces of information:

1. Type of variation (Substitution, Insertion, Silent, Intergenic, etc.) If there are multiple possibilities, annotate with the most deleterious possibility.
2. Depth of sequence coverage at the site of variation.
3. Number of reads supporting the variant.
4. Percentage of reads supporting the variant versus those supporting reference reads.
5. Allele frequency of variant from Broad Institute ExAC Project API. 
(API documentation is available here: http://exac.hms.harvard.edu/)
6. Additional optional information from ExAC that you feel might be relevant.

# Example input and output:

**INPUT:**  
**#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	normal	vaf5**  
1	100818464	.	T	C	68895.1	.	AB=0;ABP=0;AC=6;AF=1;AN=6;AO=1932;CIGAR=1X;DP=1932;DPB=1932;DPRA=0;EPP=168.744;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=70;MQMR=0;NS=2;NUMALT=1;ODDS=748.347;PAIRED=0.969979;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=77082;QR=0;RO=0;RPL=872;RPP=42.7352;RPPR=0;RPR=1060;RUN=1;SAF=1186;SAP=220.607;SAR=746;SRF=0;SRP=0;SRR=0;TYPE=snp	GT:GQ:DP:DPR:RO:QR:AO:QA	1/1/1:144.082:966:966,966:0:0:966:38541	1/1/1:144.082:966:966,966:0:0:966:38541     
1	100949860	.	G	A	49773.4	.	AB=0.468531;ABP=34.9902;AC=2;AF=0.333333;AN=6;AO=1742;CIGAR=1X;DP=3718;DPB=3718;DPRA=0;EPP=19.2103;EPPR=4.28194;GTI=1;LEN=1;MEANALT=2;MQM=70;MQMR=70;NS=2;NUMALT=1;ODDS=57.2198;PAIRED=0.947187;PAIREDR=0.966565;PAO=0;PQA=0;PQR=0;PRO=0;QA=69872;QR=77874;RO=1974;RPL=714;RPP=125.914;RPPR=360.412;RPR=1028;RUN=1;SAF=1240;SAP=681.931;SAR=502;SRF=1432;SRP=874.349;SRR=542;TYPE=snp	GT:GQ:DP:DPR:RO:QR:AO:QA	0/0/1:-0:1859:1859,871:987:38937:871:34936	0/0/1:126.955:1859:1859,871:987:38937:871:34936  
   
**OUTPUT:**  
**Chromosome	Position	Reference_Allele	Alternate_Alleles	Type_of_Variation	Depth_of_Sequence_Coverage	Allele_Count_in_Genotypes	Variant_Reads_Percentage	Allele_Frequency	Variant_Consequences	Reference_SNP_ID	Ensembl_Gene_ID	Ensembl_Transcript_ID**  
1	100818464	T	C	Single_Nucleotide_Polymorphism	1932	6	100	0.978486536522472	intron_variant,5_prime_UTR_variant	rs594529	ENSG00000079335	ENST00000455467,ENST00000542213,ENST00000336454,ENST00000361544,ENST00000370124,ENST00000370125  
1	100949860	G	A	Single_Nucleotide_Polymorphism	3718	2	33.3333	0.167615502518064	intron_variant,synonymous_variant,3_prime_UTR_variant	rs2270694	ENSG00000079335,ENSG00000228086	ENST00000542213,ENST00000432210,ENST00000336454,ENST00000361544,ENST00000544534,ENST00000370124,ENST00000370125  

# Development:

A driver R script was first developed on RStudio (under /src/annotateVCF.R). It was then built into a package (code is under /R/ExACvcf.R).

There are three steps for generate the annotated file. 

- The first step is to filter the input VCF file and get some basic information such as chromosome and position.
- The second step is to filter the INFO field and construct call commands that can be used by ExAC API bulk queries.  
- The last step is to format the array of information about the variants found in the ExAC database, and combine with the data frame generated in step one, and then return as a single data frame for output.

In order to better annotate the variants, I chose to have the following information included in the output data frame: 

Chromosome | Position | Reference_Allele | Alternate_Alleles | Type_of_Variation | Depth_of_Sequence_Coverage | Allele_Count_in_Genotypes | Variant_Reads_Percentage | Allele_Frequency | Variant_Consequences | Reference_SNP_ID | Ensembl_Gene_ID | Ensembl_Transcript_ID 


# How to test:
To test the package, you will need to do the following:

1. Download  the ExACvcf.R package. You can either download and install the tar file, or clone the GitHub repo.
   - From tar file, run: **install.packages("ExACvcf_0.0.0.9000.tar.gz")** 
   - From GitHub, run: **install_github("AlexTRee/ExACvcf")**
2. Install required packages before running ExACvcf
   - Run: **install.packages(c("httr","jsonlite"))**
3. Make sure the input file is in place. (under **/inst/extdata/input_data.vcf**)
4. Load the package, run the following commands:
   **- library(ExACvcf)**
   **- input_file <- system.file("extdata", "input_data.vcf", package = "ExACvcf")**
   **- output <- ExACvcf(input_file)**
   **- write.table(output, "./annotated_challange_data.vcf", sep = "\t", quote=FALSE)**

The output file was uploaded under (/output/annotated_challange_data.vcf)
