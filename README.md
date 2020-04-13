# ExACvcf:

ExACvcf is a prototype variant annotation tool. It annotates VCF files through Broad Institute ExAC Project and output a new VCF file annotating each variant in the file with additional information from ExAC. Each variant is annotated with the following pieces of information:

1. Type of variation (Substitution, Insertion, Silent, Intergenic, etc.) If there are multiple possibilities, annotate with the most deleterious possibility.
2. Depth of sequence coverage at the site of variation.
3. Number of reads supporting the variant.
4. Percentage of reads supporting the variant versus those supporting reference reads.
5. Allele frequency of variant from Broad Institute ExAC Project API. 
(API documentation is available here: http://exac.hms.harvard.edu/)
6. Additional optional information from ExAC that you feel might be relevant.


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
   - From tar file, run: install.packages("ExACvcf_0.0.0.9000.tar.gz") 
   - From GitHub, run: install_github("AlexTRee/TBC")
2. Install required packages before running ExACvcf
   - Run: install.packages(c("httr","jsonlite"))
3. Make sure the input file is in place. (under /inst/extdata/challenge_data.vcf)
4. Load the package, run the following commands:
   - library(ExACvcf)
   - input_file <- system.file("extdata", "challenge_data.vcf", package = "ExACvcf")
   - output <- ExACvcf(input_file)
   - write.table(output, "./annotated_challange_data.vcf", sep = "\t", quote=FALSE)

The output file was uploaded under (/output/annotated_challange_data.vcf)
