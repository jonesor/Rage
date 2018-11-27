
# regular matrices
mat_u <- rbind(c(0.1,   0,   0,   0),
               c(0.5, 0.2, 0.1,   0),
               c(  0, 0.3, 0.3, 0.1),
               c(  0,   0, 0.5, 0.6))

mat_f <- rbind(c(  0,   0, 1.1, 1.6),
               c(  0,   0, 0.8, 0.4),
               c(  0,   0,   0,   0),
               c(  0,   0,   0,   0))

mat_c <- rbind(c(  0,   0, 0.4, 0.4),
               c(  0,   0, 0.1, 0.6),
               c(  0,   0,   0,   0),
               c(  0,   0,   0,   0))


# matrices with na
mat_u_na <- rbind(c(0.1,   0,   0,   0),
                  c(0.5, 0.2, 0.1,   0),
                  c(  0, 0.3,  NA, 0.1),
                  c(  0,   0, 0.5, 0.6))

mat_f_na <- rbind(c(  0,   0, 1.1, 1.6),
                  c( NA,   0, 0.8, 0.4),
                  c(  0,   0,   0,   0),
                  c(  0,   0,   0,   0))

mat_c_na <- rbind(c(  0,   0, 0.4, 0.4),
                  c(  0,   0, 0.1, 0.6),
                  c(  0,   0,   0,   0),
                  c(  0,  NA,   0,   0))

# zero matrices
mat_u_zero <- matrix(0, nrow = 4, ncol = 4)
mat_f_zero <- matrix(0, nrow = 4, ncol = 4)
mat_c_zero <- matrix(0, nrow = 4, ncol = 4)

# singular matU
mat_u_singular <- rbind(c(0.1,   0,   0,   0),
                        c(0.5, 0.2,   0,   0),
                        c(  0, 0.3, 0.5, 0.4),
                        c(  0,   0, 0.5, 0.6))

# survivalIssue matU
mat_u_survissue <- rbind(c(0.1,   0,   0,   0),
                         c(0.5, 0.2, 0.1,   0),
                         c(  0, 0.3, 0.3, 0.1),
                         c(  0,   0, 0.7, 0.6))


# inter-reproductive stage
mat_u_inter <- rbind(c(0.1,   0,   0,   0,   0),
                     c(0.5, 0.2, 0.1,   0,   0),
                     c(  0, 0.3, 0.3, 0.1,   0),
                     c(  0,   0, 0.4, 0.4, 0.1),
                     c(  0,   0,   0, 0.1, 0.4))

mat_f_inter <- rbind(c(  0, 1.1,   0, 1.6,   0),
                     c(  0, 0.8,   0, 0.4,   0),
                     c(  0,   0,   0,   0,   0),
                     c(  0,   0,   0,   0,   0),
                     c(  0,   0,   0,   0,   0))


# age-based mpm
mat_u_age <- rbind(c(  0,   0,   0,   0),
                   c(0.8,   0,   0,   0),
                   c(  0, 0.7,   0,   0),
                   c(  0,   0, 0.6,   0))

mat_f_age <- rbind(c(  0,   0, 1.1, 1.6),
                   c(  0,   0,   0,   0),
                   c(  0,   0,   0,   0),
                   c(  0,   0,   0,   0))


# non-square
mat_u_notsq <- rbind(c(0.1,   0,   0,   0),
                     c(0.5, 0.2, 0.1,   0),
                     c(  0, 0.3, 0.3, 0.1),
                     c(  0,   0, 0.5, 0.6),
                     c(  0,   0,   0, 0.1))

mat_f_notsq <- rbind(c(  0,   0, 1.1, 1.6),
                     c(  0,   0, 0.8, 0.4),
                     c(  0,   0,   0,   0),
                     c(  0,   0,   0,   0),
                     c(  0,   0,   0,   0))


# all stages reproductive
mat_u_allrep <- rbind(c(0.1, 0.0),
                      c(0.4, 0.5))

mat_f_allrep <- rbind(c(0.7, 1.6),
                      c(  0,   0))

