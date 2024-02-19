#### Tests for small internal functions ####
test_that("Recylability testing works", {
  # expect use cases
  x <- list(a = seq(3), b = letters[seq(3)])
  expect_no_condition(
    .test_recyclable(x)
  )
  expect_true(
    .test_recyclable(x)
  )
  expect_false(
    .test_recyclable(
      list(seq(3), letters)
    )
  )

  # expect error
  x <- "list"
  expect_error(
    .test_recyclable(x),
    regexp = "`x` must be a list with vector elements"
  )
})

test_that("Recycling vectors works", {
  # expect use case
  x <- list(a = seq(3), b = letters[seq(3)])
  expect_no_condition(
    .recycle_vectors(x)
  )
  expect_identical(
    .recycle_vectors(x), x
  )
  expect_identical(
    .recycle_vectors(list(seq(3), 1)),
    list(seq(3), rep(1, 3))
  )
  # not envisaged but works to recycle lists as well
  expect_identical(
    .recycle_vectors(list(seq(3), list(1))),
    list(seq(3), as.list(rep(1, 3)))
  )
  # not envisaged as expected to be called after test for recyclability
  expect_no_condition(
    .recycle_vectors(list(seq(3), seq(2)))
  )
  expect_identical(
    .recycle_vectors(list(seq(3), seq(2))),
    list(seq(3), rep(seq(2), 3))
  )

  # expect error
  expect_error(
    .recycle_vectors("list"),
    regexp = "`x` must be a list with vector elements"
  )
  expect_error(
    .recycle_vectors(seq(3)),
    regexp = "`x` must be a list with vector elements"
  )
})
