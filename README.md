# MutationalPatterns

The MutationalPatterns R package provides a comprehensive set of flexible
functions for easy finding and plotting of mutational patterns in base
substitution catalogues.

## ChangeLog

### v0.99.6

#### Renamed functions

  * `mut_type_occurences` to `mut_type_occurrences`
  * `strand_occurences` to `strand_occurrences`

### v0.99.5

  * Added deprecation and defunct messages to functions that have
    changed since the `v0.99.0`.
  * Various small vignette and reference manual updates.

### v0.99.4

  * Internal package loading changes.
  * Removed files that do not belong to the package.

### v0.99.3

#### Renamed functions

  * `get_mut_context` to `mutation_context`
  * `get_type_context` to `type_context`
  * `get_muts` to `mutations_from_vcf`
  * `get_strand` to `strand_from_vcf`

#### Other changes

  * Added an explanation for the difference between
    SomaticSignatures and MutationalPatterns in the vignette.

### v0.99.2

#### Renamed functions

  * `vcf_to_granges` to `read_vcfs_as_granges`
  * `get_types` to `mutation_types`

### v0.99.1

#### Renamed functions
  
  * `read_vcf` to `vcf_to_granges`

#### Removed functions
  
  * `bed_to_granges`, `estimate_rank`, `rename_chrom`

#### Parameter changes
  
  * `plot_rainfall`, `vcf_to_granges`

### v0.2.0-beta

#### Updated functions
  * Faster vcf file loading
  * Automatic check and exclusion of positions in vcf with indels and/or
    multiple alternative alleles 
  * Default plotting colours to standard in the mutational signatures field

#### New features
  * Function to make chromosome names uniform according to e.g. UCSC standard
  * Transcriptional strand bias analysis
  * Signature extraction (NMF) with transcriptional strand information
  * Enrichment/depletion test for genomic annotations

Please give credit and cite MutationalPatterns R Package when you use it for
your data analysis.  A preprint of the article can be found on 
[bioRxiv](http://dx.doi.org/10.1101/071761).

# Table of Contents

* [Getting started](#getting-started)
  * [Installation](#installation)
  * [Reference genome](#reference-genome)
  * [Load data](#load-data)
  * [Make chromosome names uniform](#make-chromosome-names-uniform)
* [Mutation characteristics](#mutation-characteristics)
  * [Base substitution types](#base-substitution-types)
  * [Mutation spectrum](#mutation-spectrum)
  * [96 Mutation profile](#96-mutation-profile)
* [Mutational signatures](#mutational-signatures)
  * [De novo mutational signature extraction](#de-novo-mutational-signature-extraction-using-NMF)
  * [Fit 96 mutation profiles to known signatures](#fit-96-mutation-profiles-to-known-signatures)
* [Transcriptional strand bias](#transcriptional-strand-bias)
  * [Strand bias analysis](#strand-bias-analysis)
  * [Extract signatures with strand bias](#extract-signatures-with-strand-bias)
* [Genomic distribution](#genomic-distribution)
  * [Rainfall plot](#rainfall-plot)
  * [Enrichment or depletion of mutations in genomic regions](#enrichment-or-depletion-of-mutations-in-genomic-regions)

# Getting started

## Installation

For this package, you need R version 3.3.0 or higher.  Additionally,
you need `BiocInstaller` and `devtools` to easily load this package.

Install `BiocInstaller`:

  ```{r}
  source("https://bioconductor.org/biocLite.R")
  biocLite("BiocInstaller")
  library("BiocInstaller")
  ```

Install `devtools` to load this package:

  ```{r}
  install.packages("devtools")
  library(devtools)
  ```

Install and load `MutationalPatterns`:

  ```{r}
  options(unzip = 'internal')
  install_github("CuppenResearch/MutationalPatterns")
  library(MutationalPatterns)
  ```

## Pick a reference genome

1. List all available reference genomes (BSgenome)

  ```{r}
  library(BSgenome)
  available.genomes()
  ```
2. Download and load your reference genome of interest

  ```{r}
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
  biocLite(ref_genome)
  library(ref_genome, character.only = TRUE)
  ```
  
## Load data

Small data samples are included in the package.  You can download larger
samples from [the MutationalPatterns-data repository](https://github.com/CuppenResearch/MutationalPatterns-data/).

The data in that repository consists of somatic mutation catalogues of
nine normal human adult stem cells from three tissues (Blokzijl et al., 2016).

To load data, we need to locate it:
  ```{r}
  # Use the provides sample data
  vcf_files <- list.files(system.file("extdata", package="MutationalPatterns"),
                          pattern = ".vcf", full.names = TRUE)

  # Or list your own vcf files
  vcf_files = list.files(your_dir, pattern = ".vcf", full.names = TRUE)
  ```

And define corresponding names for the datasets:
  ```{r}
  sample_names = c("colon1", "colon2", "colon3", 
                   "intestine1", "intestine2", "intestine3",
                   "liver1", "liver2", "liver3")
  ```

This package is for the analysis of patterns in base substitution data only,
therefore indel positions and positions with multiple alternative alleles are
discarded.

Load a single VCF file:
  ```{r}
  vcf = read_vcfs_as_granges(vcf_files[1], sample_names[1], genome = "hg19")
  ```

Or load a list of VCF files:
  ```{r}
  vcfs = read_vcfs_as_granges(vcf_files, sample_names, genome = "hg19")
  ```

Include relevant metadata in your analysis, e.g. donor id, cell type, age,
tissue type, mutant or wild type
  ```{r}
  tissue = c(rep("colon", 3), rep("intestine", 3), rep("liver", 3))
  ```

## Take a subset of the chromosomes

You could, for example, select autosomal chromosomes:

  ```{r}
  auto = extractSeqlevelsByGroup(species="Homo_sapiens", 
                                 style="UCSC",
                                 group="auto")

  vcfs = lapply(vcfs, function(x) keepSeqlevels(x, auto))
  ```

# Mutation characteristics

## Base substitution types

Retrieve base substitutions from VCF object as "REF>ALT".

  ```{r}
  mutations_from_vcf(vcfs[[1]])
  ```
  
Retrieve base substitutions from vcf and convert to the 6 types of base
substitution types that are distinguished by convention: `C>A`, `C>G`, `C>T`,
`T>A`, `T>C`, `T>G`.  For example, if the reference allele is `G` and the
alternative allele is `T` (`G>T`), this functions returns the `G:C>T:A`
mutation as a `C>A` mutation.

  ```{r}
  mutation_types(vcfs[[1]])
  ```
  
Retrieve the context (1 base upstream and 1 base downstream) of the positions
in the vcf object from the reference genome.

  ```{r}
  mutation_context(vcfs[[1]], ref_genome)
  ```

Retrieve the types and context of the base substitution types for all positions
in the vcf object. For the base substitutions that are converted to the
conventional base substitution types, the reverse complement of the context is
returned.

  ```{r}
  type_context(vcfs[[1]], ref_genome)
  ```

Count mutation type occurrences for one vcf object

  ```{r}
  type_occurrences = mut_type_occurrences(vcfs[1], ref_genome)
  ```

Count mutation type occurrences for all samples in a list of vcf objects

  ```{r}
  type_occurrences = mut_type_occurrences(vcfs, ref_genome)
  ```

## Mutation spectrum

Plot mutation spectrum over all samples. Plots the mean relative contribution
of each of the 6 base substitution types. Error bars indicate standard
deviation over all samples. The n indicates the total number of mutations in
the set.

  ```{r}
  plot_spectrum(type_occurrences)
  ```

Plot mutation spectrum with distinction between C>T at CpG sites

  ```{r}
  plot_spectrum(type_occurrences, CT = TRUE)
  ```

Plot spectrum without legend

  ```{r}
  plot_spectrum(type_occurrences, CT = TRUE, legend = FALSE)
  ```

  ![spectra1](https://github.com/CuppenResearch/MutationalPatterns/blob/develop/images/spectra1.png)


Plot spectrum for each tissue separately
  ```{r}
  plot_spectrum(type_occurrences, by = tissue, CT = TRUE)
  ```

Specify 7 colors for spectrum plotting
  ```{r}
  my_colors = c("pink", "orange", "blue", "lightblue", "green", "red", "purple")
  plot_spectrum(type_occurrences, CT = T, legend = T, colors = my_colors)
  ```
  
  ![spectra2](https://github.com/CuppenResearch/MutationalPatterns/blob/develop/images/spectra2.png)

## 96 Mutation profile

Make 96 trinucleodide mutation count matrix
  ```{r}
  test_matrix = mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
  ```

Plot 96 profile of three samples
  ```{r}
  plot_96_profile(test_matrix[,c(1,4,7)])
  ```
  ![96_mutation_profile](https://github.com/CuppenResearch/MutationalPatterns/blob/develop/images/96_profile.png)

# Mutational signatures

## De novo mutational signature extraction using NMF

A critical parameter in NMF is the factorization rank, which is the number of
mutational signatures.  Determine the optimal factorization rank using the NMF
package (Gaujoux and Seoighe, 2010). As described in their paper: "...a common
way of deciding on the rank is to try different values, compute some quality
measure of the results, and choose the best value according to this quality
criteria. The most common approach is to choose the smallest rank for which
cophenetic correlation coefficient starts decreasing. Another approach is to
choose the rank for which the plot of the residual sum of squares (RSS) between
the input matrix and its estimate shows an inflection point."

  ```{r}
  # Add a tiny psuedocount to avoid a 0 in the matrix.
  mut_mat = mut_mat + 0.0001

  # Use the NMF package to generate an estimate plot.
  library("NMF")
  estimate = nmf(mut_mat, rank=2:5, method="brunet", nrun=100, seed=123456)
  plot(estimate)
  ```

  ![estim_rank](https://github.com/CuppenResearch/MutationalPatterns/blob/develop/images/estim_rank.png)

Extract and plot 2 signatures

  ```{r}
  nmf_res = extract_signatures(test_matrix, rank = 3)
  # provide signature names (optional)
  colnames(nmf_res$signatures) = c("Signature A", "Signature B")
  # plot signatures
  plot_96_profile(nmf_res$signatures)
  ```

  ![signatures](https://github.com/CuppenResearch/MutationalPatterns/blob/develop/images/signatures.png)

Plot signature contribution

  ```{r}
  # provide signature names (optional)
  rownames(nmf_res$contribution) = c("Signature A", "Signature B")
  # plot relative signature contribution
  plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "relative")
  # plot absolute signature contribution
  plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "absolute")
  ```

  ![contribution1](https://github.com/CuppenResearch/MutationalPatterns/blob/develop/images/contribution1.png)
  
  ```{r}
  # plot contribution of signatures for subset of samples with index parameter
  plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "absolute", index = c(1,2))
  # flip X and Y coordinates
  plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "absolute", coord_flip = T)
  ```
  
  ![contribution2](https://github.com/CuppenResearch/MutationalPatterns/blob/develop/images/contribution2.png)

Compare reconstructed mutation profile with original mutation profile

  ```{r}
  plot_compare_profiles(test_matrix[,1], nmf_res$reconstructed[,1], profile_names = c("Original", "Reconstructed"))
  ```

  ![originalVSreconstructed](https://github.com/CuppenResearch/MutationalPatterns/blob/develop/images/original_VS_reconstructed.png)


## Fit 96 mutation profiles to known signatures  

Download signatures from pan-cancer study (Alexandrov et al. 2013)
  
  ```{r}
  sp_url = "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
  cancer_signatures = read.table(sp_url, sep = "\t", header = T)
  # reorder (to make the order of the trinucleotide changes the same)
  cancer_signatures = cancer_signatures[order(cancer_signatures[,1]),]
  # only signatures in matrix
  cancer_signatures = as.matrix(cancer_signatures[,4:33])
  ```

Fit mutation matrix to cancer signatures. This function finds the optimal linear combination of mutation signatures that most closely reconstructs the mutation matrix by solving nonnegative least-squares constraints problem

  ```{r}
  fit_res = fit_to_signatures(test_matrix, cancer_signatures)
  # select signatures with some contribution
  select = which(rowSums(fit_res$contribution) > 0)
  # plot contribution
  plot_contribution(fit_res$contribution[select,], cancer_signatures[,select], coord_flip = F, mode = "absolute")

  ```

  ![signatures](https://github.com/CuppenResearch/MutationalPatterns/blob/develop/images/contribution_cancer_sigs.png)

Compare reconstructed mutation profile of sample 1 using cancer signatures with original profile

  ```{r}
  plot_compare_profiles(test_matrix[,1], fit_res$reconstructed[,1], profile_names = c("Original", "Reconstructed \n cancer signatures"))
  ```

  ![contribution](https://github.com/CuppenResearch/MutationalPatterns/blob/develop/images/original_VS_reconstructed_cancer_sigs.png)


# Transcriptional strand bias

## Strand bias analysis

For the mutations within genes it can be determined whether the mutation is on the transcribed or non-transcribed strand, which is interesting to study involvement of transcription-coupled repair. To this end, it is determined whether the "C" or "T" base (since by convention we regard base substitutions as C>X or T>X) are on the same strand as the gene definition. Base substitions on the same strand as the gene definitions are considered "untranscribed", and on the opposite strand of gene bodies as transcribed, since the gene definitions report the coding or sense strand, which is untranscribed. No strand information is reported for base substitution that overlap with more than one gene body.

Find gene definitions for your reference genome.

  ```{r}
  # get knowngenes table from UCSC for hg19
  biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
  library("TxDb.Hsapiens.UCSC.hg19.knownGene")
  genes_hg19 = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
  ```
  
Get transcriptional strand information for all positions in vcf 1. "Base substitions on the same strand as the gene definitions are considered. "-" for positions outside gene bodies, "U" for untranscribed/sense/coding strand, "T" for transcribed/anti-sense/non-coding strand.

  ```{r}
  strand_from_vcf(vcfs[[1]], genes_hg19)
  ```

Make mutation count matrix with transcriptional strand information (96 trinucleotides * 2 strands = 192 features). NB: only those mutations that are located within gene bodies are counted.

  ```{r}
  mut_mat_s = mut_matrix_stranded(vcfs, ref_genome, genes_hg19)
  ```
  
Perform strand bias analysis

  ```{r}
  strand_counts = strand_occurrences(mut_mat_s, by=tissue)
  strand_plot = plot_strand(strand_counts, mode = "relative")
  ```

Perform poisson test for strand asymmetry significance testing

  ```{r}
  strand_bias = strand_bias_test(strand_counts)
  strand_bias_plot = plot_strand_bias(strand_bias)
  ```
  
  ![strand](https://github.com/CuppenResearch/MutationalPatterns/blob/develop/images/strand.png)  

  
## Extract signatures with strand bias


  ```{r}
  # extract 2 signatures
  nmf_res_strand = extract_signatures(mut_mat_s, rank = 2)
  
  # provide signature names (optional)
  colnames(nmf_res_strand$signatures) = c("Signature A", "Signature B")
  # plot signatures with 192 features
  plot_192_profile(nmf_res_strand$signatures)
  
  # plot strand bias per mutation type for each signature with significance test
  plot_signature_strand_bias(nmf_res_strand$signatures)
  ```
  
  ![signatures_strand](https://github.com/CuppenResearch/MutationalPatterns/blob/develop/images/signatures_strand.png)  


# Genomic distribution

## Rainfall plot

A rainfall plot visualizes mutation types and intermutation distance. Rainfall plots can be used to visualize the distribution of mutations along the genome or a subset of chromosomes. The y-axis corresponds to the distance of a mutation with the previous mutation and is log10 transformed. Drop-downs from the plots indicate clusters or "hotspots" of mutations.

Make rainfall plot of sample 1 over all autosomal chromosomes
  ```{r}
  # define autosomal chromosomes
  chromosomes = seqnames(get(ref_genome))[1:22]
  # make rainfall plot
  plot_rainfall(vcfs[[1]], title = names(vcfs[1]), ref_genome = ref_genome, chromosomes = chromosomes, cex = 1.5)

  ```

  ![rainfall1](https://github.com/CuppenResearch/MutationalPatterns/blob/develop/images/rainfall1.png)
  
Make rainfall plot of sample 1 over chromosome 1

  ```{r}
  chromosomes = seqnames(get(ref_genome))[1]
  plot_rainfall(vcfs[[1]], title = names(vcfs[1]), ref_genome = ref_genome, chromosomes = chromosomes[1], cex = 2)
  ```
  ![rainfall2](https://github.com/CuppenResearch/MutationalPatterns/blob/develop/images/rainfall2.png)
  
  
## Enrichment or depletion of mutations in genomic regions

Test for enrichment or depletion of mutations in certain genomic regions, such as promoters, CTCF binding sites and transcription factor binding sites. To use your own genomic region definitions (based on e.g. ChipSeq experiments) specify your genomic regions in a named list of GRanges objects. Alternatively, use publically available genomic annotation data, like in the example below.

### Example: regulation annotation data from Ensembl using biomaRt

Here is an example of how to download regulation annotation data for genome build hg19. For other datasets, see biomaRt package documentation (Durinck et al. 2005). 

Install and load biomaRt package

  ```{r}
  source("https://bioconductor.org/biocLite.R")
  biocLite("biomaRt")
  library(biomaRt)
  ```

  ```{r}
  mart="ensembl"
  # list datasets available from ensembl BiomaRt
  listDatasets(useEnsembl(biomart="regulation"))
  ```

  ```{r}
  # Multicell regulatory features for hg19
  regulation_regulatory = useEnsembl(biomart="regulation", dataset="hsapiens_regulatory_feature", GRCh = 37)
  # list all possible filters
  listFilters(regulation_regulatory)
  # list all posible output attributes
  listAttributes(regulation_regulatory)
  # list all filter options for a specific attribute
  filterOptions("regulatory_feature_type_name", regulation_regulatory)
  ```

Download data from Ensembl using biomaRt

  ```{r}
  # Promoters
  promoter = getBM(attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name'), 
                   filters = "regulatory_feature_type_name", 
                   values = "Promoter", 
                   mart = regulation_regulatory)
  promoter_g = reduce(GRanges(promoter$chromosome_name, IRanges(promoter$chromosome_start, promoter$chromosome_end)))
  
    # Promoter flanking regions
    promoter_flanking = getBM(attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name'), 
                filters = "regulatory_feature_type_name", 
                values = "Promoter Flanking Region", 
                mart = regulation_regulatory)
  promoter_flanking_g = reduce(GRanges(promoter_flanking$chromosome_name, IRanges(promoter_flanking$chromosome_start, promoter_flanking$chromosome_end))) 
  
  # CTCF binding sites
  CTCF = getBM(attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name', 'cell_type_name'), 
                filters = "regulatory_feature_type_name", 
                values = "CTCF Binding Site", 
                mart = regulation_regulatory)
  # convert to GRanges object
  CTCF_g = reduce(GRanges(CTCF$chromosome_name, IRanges(CTCF$chromosome_start, CTCF$chromosome_end))) 
  
  # Open chromatin
  open = getBM(attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name'), 
                filters = "regulatory_feature_type_name", 
                values = "Open chromatin", 
                mart = regulation_regulatory)
  open_g = reduce(GRanges(open$chromosome_name, IRanges(open$chromosome_start, open$chromosome_end))) 
  

  # Transcription factor binding sites
  TF_binding = getBM(attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name'), 
                            filters = "regulatory_feature_type_name", 
                            values = "TF binding site", 
                            mart = regulation_regulatory)
  TF_binding_g = reduce(GRanges(TF_binding$chromosome_name, IRanges(TF_binding$chromosome_start, TF_binding$chromosome_end))) 
  ```

Combine all genomic regions (`GRanges` objects) in a `GRangesList`.

  ```{r}
  regions = list(promoter_g, promoter_flanking_g, CTCF_g, open_g, TF_binding_g)
  names(regions) = c("Promoter", "Promoter flanking", "CTCF", "Open chromatin", "TF binding")
  ```

## Test for significant depletion or enrichment in genomic regions

It is necessary to include a list with Granges of regions that were surveyed in your analysis for each sample, that is: positions in the genome at which you have enough high quality reads to call a mutation. This can for example be determined using CallableLoci tool by GATK. If you would not include the surveyed area in your analysis, you might for example see a depletion of mutations in a certain genomic region that is solely a result from a low coverage in that region, and therefore does not represent an actual depletion of mutations.

  ```{r}
  # Locate the file with surveyed/callable regions
  surveyed_file <- list.files(system.file("extdata",
                            package = "MutationalPatterns"),
                            pattern = ".bed",
                            full.names = TRUE)

  # Read BED file as GRanges object using rtracklayer
  library(rtracklayer)
  surveyed <- import(surveyed_file)

  # Unify to a single seqlevel style
  seqlevelsStyle(surveyed) <- "UCSC"

  # For this example we use the same surveyed file for each sample
  surveyed_list= rep(surveyed_list, 9)
  ```
  
Test for an enrichment or depletion of mutations in your defined genomic regions using a binomial test. For this test, the chance of observing a mutation is calculated as the total number of mutations, divided by the total number of surveyed bases.
  
  ```{r}
  # for each sample calculate the number of observed and expected number of mutations in each genomic regions
  distr = genomic_distribution(vcfs, surveyed_list, regions)
  # test for significant enrichment or depletion in the genomic regions
  # samples can be collapsed into groups, here the analysis is performed per tissue type
  distr_test = enrichment_depletion_test(distr, by = tissue)
  # plot enrichment depletion test results
  plot_enrichment_depletion(distr_test)
  ```

  ![distr](https://github.com/CuppenResearch/MutationalPatterns/blob/develop/images/genomic_distribution.png)
