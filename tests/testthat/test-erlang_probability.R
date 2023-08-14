test_that("Computing Erlang probabilites", {
  # using R only
  v <- prob_discrete_erlang(5, 1)

  expect_vector(v, ptype = numeric())
  expect_length(v, 12L) ## this is a known value and hardcoded

  expect_snapshot(v)
})

test_that("Computing Erlang probabilites with C++", {
  # using R only
  v <- prob_discrete_erlang_cpp(5L, 1)

  expect_vector(v, ptype = numeric())
  expect_length(v, 12L) ## this is a known value and hardcoded

  expect_snapshot(v)
})

test_that("Expect R-C++ equivalence", {
  expect_identical(
    prob_discrete_erlang_cpp(5L, 1),
    prob_discrete_erlang(5, 1)
  )
})
