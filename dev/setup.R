# setup
usethis::use_build_ignore("dev")
#usethis::use_data_raw()
devtools::document()

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
usethis::use_lgpl_license( name = "Alice Giampino and Roberto Ascari" )
usethis::use_tidy_description()

usethis::use_rcpp_armadillo()
usethis::use_package("LearnBayes")
# usethis::use_package("stats")
# usethis::use_package("spdep")
# usethis::use_pipe(export = TRUE)
# usethis::use_package("dplyr")
# usethis::use_package("sf")
# usethis::use_package("ggplot2")
# usethis::use_package("spdep")
# usethis::use_package("purrr")
# usethis::use_package("parallel")
# usethis::use_package("progress")
# usethis::use_package("gganimate")
# usethis::use_package("testthat", "Suggest")
# usethis::use_package("utils", "Suggest")
# usethis::use_package("methods")
# usethis::use_package("stats")
# usethis::use_vignette(name="epidsampler")


# Read me
usethis::use_readme_md( open = FALSE )


# Test that
# usethis::use_testthat()
# usethis::use_test("map_generation")
# usethis::use_test("map_step")


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
