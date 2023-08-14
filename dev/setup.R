# setup
usethis::use_build_ignore("dev")

# Description
usethis::use_description(
  list(
    Title = "Extended Flexible Latent Dirichlet Allocation",
    `Authors@R` = "c(
    person('Alice', 'Giampino', email = 'a.giampino@campus.unimib.it', role = c('cre', 'aut')),
    person('Roberto', 'Ascari', email = 'roberto.ascari@unimib.it', role = c('aut')))",
    Description = "Extended Flexible Latent Dirichlet Allocation.
    A tool useful for topic modelling analyses.
    Positive correlation between topic allowed.",
    URL = "https://github.com/AliceGiampino/EFLDA"
  )
)
usethis::use_lgpl_license()       # You can set another license here
usethis::use_tidy_description()   # sort fields and packages

usethis::use_readme_md( open = FALSE )

usethis::use_package("LearnBayes")
usethis::use_package("MASS")
usethis::use_package("tm")
usethis::use_package("stringr")
usethis::use_package("stats")
usethis::use_package("rlist")
usethis::use_rcpp_armadillo()
Rcpp::compileAttributes()
devtools::document()
# usethis::use_namespace(roxygen = TRUE)

# Tests
# usethis::use_testthat()
# usethis::use_test("gd")


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

# Vignette
# usethis::use_vignette(name="linear_model_with_gradient_package")
# Install your package as developer: devtools::install(build_vignettes = TRUE)
# Install your package from github: remotes::install_github('andreamelloncelli/myShapes', force = T, build_vignettes = T)
# Install your pagkage from file install.packages('.', repos = NULL, force = T, build_vignettes = T)

# Deploy
# devtools::missing_s3()
# devtools::check()
# rhub::check_for_cran()








