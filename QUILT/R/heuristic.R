make_rhb_t_local <- function(
    which_haps_to_use,
    hapMatcherR,
    hapMatcher,
    use_hapMatcherR,
    distinctHapsB,
    eMatDH_special_matrix,
    eMatDH_special_matrix_helper,
    nGrids
) {
    ## build eHapsCurrent_tc here
    Ksubset <- length(which_haps_to_use)
    rhb_t_local <- array(0L, c(Ksubset, nGrids))
    for(iGrid in 1:nGrids) {
        for(i_k in 1:Ksubset) {
            k <- which_haps_to_use[i_k]
            if (use_hapMatcherR) {
                i <- as.integer(hapMatcherR[k, iGrid])
            } else {
                i <- as.integer(hapMatcher[k, iGrid])
            }
            if (i > 0) {
                ## note, could make distinctHapsB pretty easily if needed for RAM
                rhb_t_local[i_k, iGrid] <- distinctHapsB[i, iGrid]
            } else {
                rhb_t_local[i_k, iGrid] <- rcpp_simple_binary_matrix_search(
                    val = k - 1,
                    mat = eMatDH_special_matrix,
                    s1 = eMatDH_special_matrix_helper[iGrid, 1],
                    e1 = eMatDH_special_matrix_helper[iGrid, 2]
                )
            }
        }
    }
    rhb_t_local
}
    


compare_heuristic_approaches <- function(
    hapProbs_t,
    which_haps_to_use_quilt1,
    which_haps_to_use_zilong_A,
    which_haps_to_use_zilong_B,
    which_haps_to_use_mspbwt_A,
    which_haps_to_use_mspbwt_B,
    hapMatcherR,
    hapMatcher,
    use_hapMatcherR,
    distinctHapsB,
    eMatDH_special_matrix,
    eMatDH_special_matrix_helper,
    outputdir,
    sample_name,
    regionName,
    i_gibbs_sample,
    i_it,
    nGrids,
    mspbwtL,
    mspbwtM
) {

    filename <- file.path(
        outputdir,
        "plots",
        paste0("heuristics.", sample_name, ".",regionName, ".", i_gibbs_sample, ".", i_it, ".pdf")
    )

    ## save(
    ## hapProbs_t,
    ## which_haps_to_use_zilong_A,
    ## which_haps_to_use_zilong_B,
    ## which_haps_to_use_quilt1,
    ## which_haps_to_use_mspbwt,
    ## use_hapMatcherR,
    ## distinctHapsB,
    ## eMatDH_special_matrix,
    ## eMatDH_special_matrix_helper,
    ## outputdir,
    ## sample_name,
    ## regionName,
    ## i_gibbs_sample,
    ## i_it,
    ## nGrids,
    ## file = paste0(filename, ".RData"))
    ## print(paste0("SAVING FOR FILE:", paste0(filename, ".RData")))
    
    pdf(filename, height = 30, width = 8)
    par(mfrow = c(5, 2))
    
    ## so suppose we had rhb_t
    ## then would want for each of these a row against current hapProbs_t
    ## showing agree vs disagree
    for(i in 1:5) {
        
        if (i == 1) {
            which_haps_to_use <- which_haps_to_use_quilt1
            approach <- "QUILT1"
        } else if (i == 2) {
            which_haps_to_use <- which_haps_to_use_zilong_A
            approach <- "Zilong A"
        } else if (i == 3) {
            which_haps_to_use <- which_haps_to_use_zilong_B
            approach <- "Zilong B"
        } else if (i == 4) {
            which_haps_to_use <- which_haps_to_use_mspbwt_A
            approach <- "mspbwt A"
        } else if (i == 5) {
            which_haps_to_use <- which_haps_to_use_mspbwt_B
            approach <- "mspbwt B"
        }


        Ksubset <- length(which_haps_to_use)

        if (Ksubset > 0) {
            
            rhb_t_local <- make_rhb_t_local(
                which_haps_to_use,
                hapMatcherR,
                hapMatcher,
                use_hapMatcherR,
                distinctHapsB,
                eMatDH_special_matrix,
                eMatDH_special_matrix_helper,
                nGrids
            )
            
            for(i_hap in 1:2) {

                hap <- rcpp_int_contract(round(hapProbs_t[i_hap, ]))
                what_is_plotted <- array(FALSE, dim(rhb_t_local))
                
                plot(x = 0, y = 0, xlim = c(0, nGrids + 1), ylim = c(1, Ksubset), axes = TRUE, xlab = "Grid", ylab = "Hap match", main = paste0(approach, ", hap ", i_hap))
                xleft <- 1:nGrids - 0.5
                xright <- xleft + 1
                
                for(i_k in 1:Ksubset) {

                    ## plot upside down
                    x <- rhb_t_local[Ksubset - i_k + 1, ] == hap
                    ## censor if fewer than 5 in a row?
                    x2 <- rep(FALSE, nGrids)
                    a <- rle(x)
                    keep <- which((a$values == TRUE) & (a$lengths > 10))
                    b <- cumsum(a$lengths)
                    starts <- c(1, b[-length(b)] + 1)
                    ends <- b
                    if (length(keep) > 0) {
                        for(i in keep) {
                            x2[starts[i]:ends[i]] <- TRUE
                        }
                    }
                    what_is_plotted[i_k, ] <- x2

                    col <- rep("white", nGrids)
                    col[x2] <- "black"
                    rect(
                        xleft = xleft,
                        xright = xright,
                        ybottom = rep(i_k, nGrids) - 0.5,
                        ytop = rep(i_k, nGrids) + 0.5,
                        col = col,
                        border = NA
                    )
                }
            }

        }

    }

    dev.off()

    
}
