
# Setup -------------------------------------------------------------------

# avoid problem with the dev scripts: dev/package-utility.R (this file)
dir.create("dev")
# save this file in `dev` as `setup.R`
usethis::use_build_ignore("dev")
# now you can save or move this file in "dev"

# Description
usethis::use_description(
  list(
    Title = "Gradient",
    `Authors@R` = "c(
    person('Alice', 'Giampino', email = 'a.giampino@campus.unimib.it', role = c('cre', 'aut')),
    person('Roberto', 'Ascari', email = 'roberto.ascari@unimib.it', role = c('aut')))",
    Description = "Extended Flexible Latent Dirichlet Allocation",
    URL = "https://github.com/AliceGiampino/EFLDA"
  )
)
usethis::use_lgpl_license()       # You can set another license here
usethis::use_tidy_description()   # sort fields and packages

# Read me
usethis::use_readme_md( open = FALSE )
# usethis::use_code_of_conduct()
# usethis::use_lifecycle_badge( "Experimental" )
# usethis::use_news_md( open = FALSE )


## Use tests: if you want to use tests
# usethis::use_testthat()
## Reinstall devtools from R 3.x to R 4.x
# install.packages(c("devtools", "pkgload"), force = T)

# Develop -----------------------------------------------------------------

usethis::use_git_config(
  scope = "user",
  user.name = "AliceGiampino",
  user.email = "giampinoalice@gmail.com"
)
usethis::use_git()


# CI
# usethis::use_travis()
# usethis::use_coverage()
# usethis::use_lifecycle_badge("experimental")
# usethis::use_github_action_check_standard(save_as="test.yaml")
# usethis::use_github_action("test-coverage")

# Website
# usethis::use_pkgdown()
# pkgdown::build_site()
# usethis::use_github_action("pkgdown")

## Add a package
usethis::use_rcpp_armadillo()
usethis::use_package("LearnBayes")
usethis::use_package("MASS")
# usethis::use_package( "dplyr" )
# remeber to add it to ROXYGEN or NAMESPACE:
#' @import dplyr  # ROXYGEN
#' import(dplyr)  # NAMESPACE

## If you want to use roxygen, enable ROXYGEN in the project.
# Menu: tools > Project options > build tools > generate the documentation with roxygen
usethis::use_namespace(roxygen = TRUE)
devtools::document() # to fill NAMESPACE and documentation with ROXYGEN comments
# or roxygen2::roxygenise() # converts roxygen comments to .Rd files.
# or [Ctrl + Shift + D] in RStudio

# insert the documentation over a function definition Ctrl + Shift + Alt + r


## Build or load
# Load the package [CTRL + SHIFT + L] or install-and-reload [CTRL + SHIFT + B]

## Check the package for Cran or [CTRL + SHIFT + E]
devtools::check(document = FALSE) # check the package

## Add internal datasets
## If you want to provide data along with your package
# usethis::use_data_raw( name = "my_dataset", open = FALSE )

## Tests
## Add one line by test you want to create
# usethis::use_test( "division" )
# Test that
# usethis::use_testthat()
# usethis::use_test("map_generation")
# usethis::use_test("map_step")

## Vignette
# usethis::use_vignette("ThisTidyPackage")
# devtools::build_vignettes()
# Install the package and see it with `vignette("ThisTidyPackage")`
# List your vignettes: vignette(package = 'myShapes')
# Install your package as developer: devtools::install(build_vignettes = TRUE)
# Install your package from github: remotes::install_github('andreamelloncelli/myShapes', force = T, build_vignettes = T)
# Install your pagkage from file install.packages('.', repos = NULL, force = T, build_vignettes = T)

# Deploy ------------------------------------------------------------------

# devtools::missing_s3()
#
# devtools::check()
# rhub::check_for_cran()
