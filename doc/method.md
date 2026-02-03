The CeleScope-atac(v1.8.1; [https://github.com/singleron-RD/CeleScope_ATAC](https://github.com/singleron-RD/CeleScope_ATAC)) analysis pipeline was employed for initial processing of scATAC-seq data.  Raw sequencing reads underwent demultiplexing, alignment, duplicate marking, peak calling, and cell calling to produce the peak matrix. Ensembl release 109 human and Ensembl release 109 mus musculus were used as reference genome. Alignment were implemented through Chromap(v0.2.6; https://github.com/haowenz/chromap/blob/master/README.md). Peak-calling were implemented throught Macs2(v2.2.7.1; https://github.com/macs3-project/MACS/blob/master/README.md)

## Reference

- Zhang, H., Song, L., Wang, X. et al. Fast alignment and preprocessing of chromatin profiles with Chromap. Nat Commun 12, 6566 (2021). https://doi.org/10.1038/s41467-021-26865-w

- Zhang, Y., Liu, T., Meyer, C.A. et al. Model-based Analysis of ChIP-Seq (MACS). Genome Biol 9, R137 (2008). https://doi.org/10.1186/gb-2008-9-9-r137
