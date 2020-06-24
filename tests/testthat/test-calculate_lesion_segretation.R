context("test-calculate_lesion_segregation")

# To test mut_matrix, we need to load the reference genome first.
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

# Load GRangesList
grl <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
                             package="MutationalPatterns"))
sample_names <- c(
   "colon1", "colon2", "colon3",
   "intestine1", "intestine2", "intestine3",
   "liver1", "liver2", "liver3")


#Perform lesion segregation calculations
output = calculate_lesion_segregation(grl, sample_names)
output_per_type = calculate_lesion_segregation(grl, sample_names, 
                                               split_by_type = TRUE, ref_genome = ref_genome)
output_walf = calculate_lesion_segregation(grl,
                                           sample_names,
                                           test = "walf-wolfowitz")

test_that("Output has correct class",{
    expect_true(inherits(output, c("tbl_df")))
    expect_true(inherits(output_per_type, c("tbl_df")))
    expect_true(inherits(output_walf, c("tbl_df")))
    
})

test_that("Output has correct dimensions",{
    expect_equal(dim(output), c(9, 8))
    expect_equal(dim(output_per_type), c(9,8))
    expect_equal(dim(output_walf), c(9,5))
})

expected <- readRDS(system.file("states/lesion_segregation.rds",
                                package="MutationalPatterns"))
test_that("transforms correctly", {
    expect_equal(output, expected)
})
