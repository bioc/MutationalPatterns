##
## Removed functions
##

bed_to_granges <- function(bed_files, region_names)
{
    .Defunct("import", package="rtracklayer",
            msg = paste("'bed_to_granges' has been removed.\nUse 'import'",
                        "from package 'rtracklayer' instead.\n\nThe",
                        "following snippet replaces this function:\n",
                        "\n    library(rtracklayer)",
                        "\n    bed <- import(<your file>)",
                        "\n    seqlevelsStyle(bed) <- \"UCSC\"\n"))
}

estimate_rank <- function(mut_matrix, rank_range, nrun=100)
{
    .Defunct("nmf", package="NMF",
            msg = paste("'estimate_rank' has been removed.\nUse 'nmf'",
                        "from package 'NMF' instead.\n\nThe following",
                        "snippet replaces this function:\n",
                        "\n    library(NMF)",
                        "\n    estimate <- nmf(<mut_matrix>,",
                        "rank=<rank_range>, method=\"brunet\", nrun=<nrun>,",
                        "seed=123456)",
                        "\n    plot(estimate)\n"))
}
