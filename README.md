Pipeline for the analysis of the linkage disequilibrium between deleterious in a population with 182 *Capsella grandiflora* individuals (data from Josephs et al., 2015) and in 33 *Capsella orientalis* individuals (data from Agren et al., 2014; Huang et al., 2018; Kryvokhyzha et al., 2019). The list of accessions can be found in `Accession_SRA_Number/ID_ind`.

To run the pipeline correctly, shell scripts should be executed in the order indicated by their number. 

To run the scripts in `3_Mapping/` the file `picard.jar` (available [here](https://github.com/broadinstitute/picard/releases/tag/3.0.0)) needs to be added in `Picard/`.

To run the scripts in `3_Mapping/`, `4_Variant_Calling/` and `5_Building_Sift_DB/` the reference genome Cr145 of *Capsella rubella* `Cr145.fasta` needs to be added to `Reference_genome/`.

To run the scripts in `5_Building_Sift_DB/` the annotation of the reference genome of *Capsella rubella* `Cr145.gff`

References:

