test_that("vital_rates functions work correctly", {
  vrOut <- vital_rates(matU = mat_u, matF = mat_f)
  
  expect_type(vrOut, "list")
  expect_length(vrOut, 5)
  expect_equal(vrOut$clo, 0)
  
  vrOut_c <- vital_rates(matU = mat_u,
                        matF = mat_f,
                        matC = mat_c)
  
  expect_type(vrOut_c, "list")
  expect_length(vrOut_c, 5)
  
  
  vrOut_ss <-
    vital_rates(matU = mat_u,
               matF = mat_f,
               splitStages = "all")
  
  expect_type(vrOut_ss, "list")
  expect_length(vrOut_ss, 5)
  expect_equal(vrOut_ss$clo, 0)
  
  vrOut_ss2 <-
    vital_rates(matU = mat_u,
               matF = mat_f,
               splitStages = "ontogeny")
  
  expect_type(vrOut_ss2, "list")
  expect_length(vrOut_ss2, 9)
  expect_equal(vrOut_ss2$cloAdu, 0)
  
  vrOut_ss3 <-
    vital_rates(
      matU = mat_u,
      matF = mat_f,
      splitStages = "matrixStages",
      matrixStages = c("active", "active", "active", "dorm")
    )
  
  expect_type(vrOut_ss3, "list")
  expect_length(vrOut_ss3, 10)
  expect_equal(vrOut_ss3$clo, 0)
  
  vrOut_wt <- vital_rates(matU = mat_u,
                         matF = mat_f,
                         weights = "SSD")
  
  expect_type(vrOut_wt, "list")
  expect_length(vrOut_wt, 5)
  expect_equal(vrOut$clo, 0)
  
  vrOut_wt2 <-
    vital_rates(matU = mat_u,
               matF = mat_f,
               weights = c(1, 2, 3, 4))
  
  expect_type(vrOut_wt2, "list")
  expect_length(vrOut_wt2, 5)
  expect_equal(vrOut_wt2$clo, 0)
})


test_that("vital_rates functions warn and fail gracefully", {
  expect_error(vital_rates(mat_u_na, mat_f_na))
  
  expect_error(vital_rates(mat_u, mat_f_na))
  
  expect_error(vital_rates(mat_u_na, mat_f))
  
  expect_error(vital_rates(mat_u, mat_f, weights = 3))
  
  expect_error(vital_rates(mat_u, mat_f, splitStages =  "wrongName"))
  
  expect_error(vital_rates(
    mat_u,
    mat_f,
    splitStages = "matrixStages",
    matrixStages = c("active", "active", "active", "x")
  ))
  
  expect_error(vital_rates(mat_u, mat_f, splitStages = "matrixStages"))
  
  expect_error(vital_rates(
    mat_u,
    mat_f,
    splitStages = "matrixStages",
    matrixStages = c("active", "active", "active")
  ))
  
  
})
