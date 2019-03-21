# what are the raw data?  
16s based reads (in .fq format), with forward primers   
(fecal samples -> microbial genomic DNA -> PCR (338F, 806R) -> ribosomal V3&V4 rRNA genes -> Miseq PE300)  

# how are the data processed?  
R packages involved: dada2 (https://www.bioconductor.org/packages/release/bioc/manuals/dada2/man/dada2.pdf), phyloseq and ggplot2  

# which pipeline of dada2? 
DADA2 ITS workflow: http://benjjneb.github.io/dada2/ITS_workflow.html  

# credits  
Dr. Xinyue Zhang, Dr. Ting Zhang @ Shanghai Children's Hospital: for providing IBD samples  

