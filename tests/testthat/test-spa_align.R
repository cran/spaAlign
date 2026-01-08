test_that("spa_align preserves reference depth", {
  ref <- data.frame(Depth_m = 1:10, A = rnorm(10))
  plugs <- data.frame(Depth_m = c(2, 5, 7),
                      TOC = c(3.0, 3.3, 3.5))

  out <- spa_align(
    ref,
    plugs = plugs,
    depth_col = "Depth_m",
    rule = 1,
    add_suffix = TRUE
  )

  expect_equal(nrow(out), nrow(ref))
  expect_equal(max(out$Depth_m), max(ref$Depth_m))
  expect_true(all(is.na(out$TOC_plugs[8:10])))
})

test_that("spa_align handles duplicate depths deterministically", {
  ref <- data.frame(Depth_m = 1:6)
  plugs <- data.frame(Depth_m = c(2, 2, 4, 5),
                      TOC = c(3.0, 3.1, 3.4, 3.6))

  out <- spa_align(ref, plugs = plugs, depth_col = "Depth_m", rule = 1, add_suffix = TRUE)

  expect_equal(nrow(out), nrow(ref))
  expect_true(all(is.finite(out$Depth_m)))
})
