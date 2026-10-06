# create the release issue to keep track of the tasks needed for submission
use_release_issue()

# first release needs a comment for CRAN
use_cran_comments()

#
check(remote = TRUE)


check_win_devel()
rhub::rhub_check()
submit_cran()
