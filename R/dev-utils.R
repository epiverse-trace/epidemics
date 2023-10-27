# This unexported function adds a custom item to `usethis::use_release_issue()`
release_bullets <- function() {
  c(
    "Run `goodpractice::gp()`",
    "Review [WORDLIST](https://docs.cran.dev/spelling#wordlist)",
    "Check if `# nolint` comments are still needed with recent lintr releases"
  )
}

# lintr does not recognise the following variables used in a non-standard way
# with cli, so they have to be declared as global variables here.
# This fix is suggested in the following issue:
# https://github.com/r-lib/lintr/issues/358#issuecomment-535991461
globalVariables(
  c(
    "header",
    "name"
  )
)
