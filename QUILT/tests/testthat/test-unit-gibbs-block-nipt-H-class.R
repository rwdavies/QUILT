test_that("can sample H given H_class in Rcpp", {

    set.seed(101)
    
    ff <- 0.2
    ## check approx frequencies between them
    H_class <- rep(as.integer(0:7), each = 1e5) ## a lot!
    ## 
    H_R <- sample_H_using_H_class(H_class, ff)
    H_Rcpp <- integer(length(H_class))
    rcpp_sample_H_using_H_class(H_class, H_Rcpp, ff)

    ## check, expect normal p-value for t.test
    for(i in 0:7) {
        w <- H_class == i
        a <- H_R[w]
        x <- c(sum(a == 1), sum(a == 2), sum(a == 3))
        a <- H_Rcpp[w]
        y <- c(sum(a == 1), sum(a == 2), sum(a == 3))
        ## print(rbind(x, y))
        p_value <- t.test(rbind(x, y))$p.value
        expect_true(p_value > 1e-4)
    }
        
    ## print(table(H_class, H_R))
    ## print(table(H_class, H_Rcpp))
    

})


test_that("can use ns for H_class checking", {

    rr <- rbind(
        c(1, 2, 3),
        c(1, 3, 2),
        c(2, 1, 3),
        c(2, 3, 1),
        c(3, 1, 2),
        c(3, 2, 1)
    )
    rx <- rbind(
        c(1, 2, 3),
        c(1, 3, 2),
        c(2, 1, 3),
        c(3, 1, 2),
        c(2, 3, 1),
        c(3, 2, 1)
    )
    ## about how to apply
    ## e.g. ir_chosen = 4
    ## is rr 2, 3, 1
    ## meaning
    ## the new 1 is the old 2
    ## the new 2 is the old 3
    ## the new 3 is the old 1
    ## so for instance for H_old
    ## H_new <- c(3, 1, 2)[H_old]
    ## so e.g. the new count of 1s is the same as the old count of 2s
    
    set.seed(9119)
    ## do in raw form
    ## probabilities are not relevant I believe
    vals <- cbind(
        sample(c(0, 1), 501, replace = TRUE),
        sample(c(0, 1), 501, replace = TRUE),
        sample(c(0, 1), 501, replace = TRUE)
    )
    ## 0 = 0 0 0
    ## 1 = 1 0 0
    ## 2 = 0 1 0
    ## 3 = 0 0 1
    ## 4 = 1 1 0
    ## 5 = 1 0 1
    ## 6 = 0 1 1
    ## 7 = 1 1 1
    get_H_class <- function(vals) {
        apply(vals, 1, function(x) {
            if (sum(x == 0) == 3) {
                H_class <- 0
            } else if (sum(x == 0) == 2) {
                H_class <- which.max(x == 1)
            } else if (sum(x == 0) == 1) {
                H_class <- 7 - which.max(x == 0)
            } else {
                H_class <- 7
            }
            H_class
        })
    }
    H_class <- get_H_class(vals)
    n <- c(
        sum(H_class == 1),
        sum(H_class == 2),
        sum(H_class == 3),
        sum(H_class == 4),
        sum(H_class == 5),
        sum(H_class == 6)
    )

    for(ir_chosen in 1:6) {

        cur_H_class <- get_H_class(vals[, rr[ir_chosen, ]])
        
        n_alt1 <- c(
            sum(cur_H_class == 1),
            sum(cur_H_class == 2),
            sum(cur_H_class == 3),
            sum(cur_H_class == 4),
            sum(cur_H_class == 5),
            sum(cur_H_class == 6)
        )

        ## check this is right manually
        expect_equal(n[1:3][rr[ir_chosen, ]], n_alt1[1:3])

        ## note, this code is used later, so we are validating the code rather than function
        n_alt2 <- c(
            n[rr[ir_chosen, 1]],
            n[rr[ir_chosen, 2]],
            n[rr[ir_chosen, 3]],
            n[7 - rr[ir_chosen, 3]],
            n[7 - rr[ir_chosen, 2]],
            n[7 - rr[ir_chosen, 1]]
        )

        expect_equal(n_alt1, n_alt2)

        ## also make sure I can do this manually
        one_based_swap <- c(1, 1 + rx[ir_chosen, 1:3], 8 - rx[ir_chosen, 3:1], 8)        
        H_class_new <- one_based_swap[H_class + 1] - 1

        n_alt3 <- c(
            sum(H_class_new == 1),
            sum(H_class_new == 2),
            sum(H_class_new == 3),
            sum(H_class_new == 4),
            sum(H_class_new == 5),
            sum(H_class_new == 6)
        )
        
        expect_equal(n_alt1, n_alt3)
        
    }

    ## AWW YEAH BABY looks good

})

    ## sapply(1:1000, function(i) {
    ##     x <- sample(c("none", "one", "three"), 1)
    ##     y <- c(0, 0, 0)
    ##     if (x == "one") {
    ##         y[sample(c(1, 2, 3), prob = c(0.5, 0.5 - ff / 2, ff / 2))] <- 1
    ##     } else if (x == "one") {
    ##         ## ONE of them is unique, choose is
    ##         l1 <- sample(1:3)
    ##         if (l1 == 1) { ## unique is mat trans
    ##             l <- sample(c(1,
    ##         }
    ##         l <- sample(c(1, 2, 3), prob = c(0.5, 0.5 - ff / 2, ff / 2))
    ##     }
    ## })

