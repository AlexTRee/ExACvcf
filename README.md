# Challenge description:

For this challenge, you are asked to prototype a variant annotation tool. We will provide you with a VCF file, and you will create a small software program to output a table annotating each variant in the file. Each variant must be annotated with the following pieces of information:

1. Type of variation (Substitution, Insertion, Silent, Intergenic, etc.) If there are multiple possibilities, annotate with the most deleterious possibility.
2. Depth of sequence coverage at the site of variation.
3. Number of reads supporting the variant.
4. Percentage of reads supporting the variant versus those supporting reference reads.
5. Allele frequency of variant from Broad Institute ExAC Project API. 
(API documentation is available here: http://exac.hms.harvard.edu/)
6. Additional optional information from ExAC that you feel might be relevant.

For this project please upload all relevant code (written in whatever language you like) along with the annotated VCF file to a Github account and provide the link to the below email address. Please note that work will be assessed based on quality of code and documentation more-so than the annotation.

