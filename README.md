Pipeline for the analysis of the linkage disequilibrium between deleterious in a population with 182 *Capsella grandiflora* individuals (data from Josephs et al., 2015) and in 33 *Capsella orientalis* individuals (data from Agren et al., 2014; Huang et al., 2018; Kryvokhyzha et al., 2019). The list of accessions can be found in `Accession_SRA_Number/ID_ind`.

To run the pipeline correctly, shell scripts should be executed in the order indicated by their number. 

To run the scripts in `3_Mapping/` the file `picard.jar` (available [here](https://github.com/broadinstitute/picard/releases/tag/3.0.0)) needs to be added in `Picard/`.

To run the scripts in `3_Mapping/`, `4_Variant_Calling/` and `5_Building_Sift_DB/` the reference genome Cr145 of *Capsella rubella* `Cr145.fasta` needs to be added to `Reference_genome/`.

To run the scripts in `5_Building_Sift_DB/` the annotation of the reference genome of *Capsella rubella* `Cr145.gff` needs to be added to `Reference_genome/` and the protein database `uniref100.fasta` (available [here](https://www.uniprot.org/help/downloads)) to `Protein_DB/`.

To run the scripts in `6_LD_Analysis` the sequence of *Neslia paniculata* needs to be added to `Neslia_sequence/`.

References:

<div class="csl-entry">Ågren, J. A., Wang, W., Koenig, D., Neuffer, B., Weigel, D., &#38; Wright, S. I. (2014). Mating system shifts and transposable element evolution in the plant genus Capsella. <i>BMC Genomics</i>, <i>15</i>(1), 602. https://doi.org/10.1186/1471-2164-15-602</div>

<div class="csl-entry">Huang, H. R., Liu, J. J., Xu, Y., Lascoux, M., Ge, X. J., &#38; Wright, S. I. (2018). Homeologue-specific expression divergence in the recently formed tetraploid Capsella bursa-pastoris (Brassicaceae). <i>New Phytologist</i>, <i>220</i>(2), 624–635. https://doi.org/10.1111/nph.15299</div>

<div class="csl-entry">Josephs, E. B., Lee, Y. W., Stinchcombe, J. R., &#38; Wright, S. I. (2015). Association mapping reveals the role of purifying selection in the maintenance of genomic variation in gene expression. <i>Proceedings of the National Academy of Sciences of the United States of America</i>, <i>112</i>(50), 15390–15395. https://doi.org/10.1073/pnas.1503027112</div>

<div class="csl-entry">Kryvokhyzha, D., Milesi, P., Duan, T., Orsucci, M., Wright, S. I., Glémin, S., &#38; Lascoux, M. (2019). Towards the new normal: transcriptomic convergence and genomic legacy of the two subgenomes of an allopolyploid weed (Capsella bursa-pastoris). <i>PLoS Genetics</i>, <i>15</i>(5), e1008131. https://doi.org/10.1371/journal.pgen.1008131</div>
