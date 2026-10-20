# create the release issue to keep track of the tasks needed for submission
usethis::use_release_issue()

# first release needs a comment for CRAN
usethis::use_cran_comments()

# checks the package as "CRAN incoming", e.g. among other things it checks
# the URLs (not usually checked when running devtools::check())
devtools::check(remote = TRUE)

# checks the package with win-builder
devtools::check_win_devel()

# rhub GHA workflow is setup, if it isn't then set it up with rhub::rhub_setup()
# check the rhub setup is correct with rhub::rhub_doctor()
rhub::rhub_check() # I chose 1-5 platforms to run the check on

# submit to CRAN (an email will be sent to validate submission)
devtools::submit_cran()
