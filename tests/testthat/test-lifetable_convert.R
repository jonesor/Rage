
test_that("surv_conversion functions work correctly", {

  # convert from lx
  lx1 <- c(1, 0.8, 0.7, 0.5, 0.3, 0.1)
  px1 <- lx_to_px(lx1)
  hx1 <- lx_to_hx(lx1)

  # convert from px
  lx2 <- px_to_lx(px1)
  hx2 <- px_to_hx(px1)

  expect_true(all(lx2 == lx1))
  expect_true(all(hx2 == hx1 | (is.na(hx2) & is.na(hx1))))

  # convert from hx
  lx3 <- hx_to_lx(hx2)
  px3 <- hx_to_px(hx2)

  expect_true(all(lx3 == lx1))
  expect_true(all(px3 == px1 | (is.na(px3) & is.na(px3))))
})


test_that("surv_conversion functions warn and fail gracefully", {
  expect_true(is.na(lx_to_px(1)))
  expect_true(is.na(lx_to_hx(1)))
})
