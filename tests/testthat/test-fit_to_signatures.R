context("test-fit_to_signatures")

#Get mut_mat
mut_mat = readRDS(system.file("states/mut_mat_data.rds",
                              package="MutationalPatterns"))

#Get signatures
filename <- system.file("extdata/snv_signatures_probabilities.txt",
                        package="MutationalPatterns")
signatures <- read.table(filename, sep = "\t", header = TRUE)
signatures = as.matrix(signatures[,-c(1,2)])


output = fit_to_signatures(mut_mat, signatures)

expected <- readRDS(system.file("states/snv_refit.rds",
                                package="MutationalPatterns"))

test_that("Output has correct class",{
    expect_true(inherits(output, "list"))
    expect_true(inherits(output$contribution, "matrix"))
    expect_true(inherits(output$reconstructed, "matrix"))

})

test_that("Output is equal to expected", {
    expect_equal(output, expected)
})

#Get indel mut_mat
indel_counts = readRDS(system.file("states/blood_indel_counts.rds", package = "MutationalPatterns"))
indel_m = indel_counts %>% 
    dplyr::select(-muttype, -muttype_sub) %>% 
    as.matrix()

#Get indel signatures
filename <- system.file("extdata/indel_signatures_probabilities.txt",
                        package="MutationalPatterns")
signatures <- read.table(filename, sep = "\t", header = TRUE)
signatures = as.matrix(signatures[,-c(1)])

expected <- readRDS(system.file("states/indel_refit.rds",
                                package="MutationalPatterns"))

test_that("Refitting indels gives expected output.", {
    output = fit_to_signatures(indel_m, signatures)
    expect_equal(output, expected)
})

#Get dbs mut_mat
dbs_counts = readRDS(system.file("states/blood_dbs_counts.rds", package = "MutationalPatterns"))
dbs_m = dbs_counts %>% 
    dplyr::select(-REF, -ALT) %>% 
    as.matrix()

#Get dbs signatures
filename <- system.file("extdata/dbs_signatures_probabilities.txt",
                        package="MutationalPatterns")
signatures <- read.table(filename, sep = "\t", header = TRUE)
signatures = as.matrix(signatures[,-c(1)])

expected <- readRDS(system.file("states/dbs_refit.rds",
                                package="MutationalPatterns"))

test_that("Refitting dbss gives expected output.", {
    output = fit_to_signatures(dbs_m, signatures)
    expect_equal(output, expected)
})
