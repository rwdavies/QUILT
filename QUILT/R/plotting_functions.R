plot_1_dosage_vs_truth <- function(dosage, truth, ancAlleleFreqAll, inRegion2, smoothV, Ls, scale, ybottom, col, minAFForBounding = 0.01) {
    ##
    ancAlleleFreqAll[ancAlleleFreqAll < minAFForBounding] <- minAFForBounding
    ancAlleleFreqAll[ancAlleleFreqAll > (1 - minAFForBounding)] <- (1 - minAFForBounding)
    ## 
    dosage2 <- (dosage - 1 * ancAlleleFreqAll) / sqrt(1 * ancAlleleFreqAll * (1 - ancAlleleFreqAll))
    truth[is.na(truth)] <- ancAlleleFreqAll[is.na(truth)]
    truth2 <- (truth - 1 * ancAlleleFreqAll) / sqrt(1 * ancAlleleFreqAll * (1 - ancAlleleFreqAll))
    ## difference!
    diff <- smooth_vector(x = abs(dosage2 - truth2), n = smoothV)
    diff <- diff / 3
    diff[diff > 1] <- 1
    r2 <- round(cor(dosage2[inRegion2], truth2[inRegion2], method = "pearson", use = "pairwise.complete") ** 2, 3)    
    print(r2)
    points(x = Ls, y = ybottom + scale * diff, col = col, type = "l")
    text(x = Ls[1], y = ybottom + scale * 1, labels = paste0("r2 = ", r2), pos = 1, cex = 1.5)
    return(r2)
}
plot_2_dosage_vs_truth <- function(dosage, truth, ancAlleleFreqAll, inRegion2, smoothV, Ls, scale, ybottom, col, minAFForBounding = 0.01) {
    ##
    ancAlleleFreqAll[ancAlleleFreqAll < minAFForBounding] <- minAFForBounding
    ancAlleleFreqAll[ancAlleleFreqAll > (1 - minAFForBounding)] <- (1 - minAFForBounding)
    ## 
    dosage2 <- dosage - 2 * ancAlleleFreqAll
    truth[is.na(truth)] <- ancAlleleFreqAll[is.na(truth)]
    truth2 <- truth - 2 * ancAlleleFreqAll
    ## difference!
    diff <- smooth_vector(x = abs(dosage2 - truth2), n = smoothV)
    diff <- diff / 3
    diff[diff > 1] <- 1
    r2 <- round(cor(dosage2[inRegion2], truth2[inRegion2], method = "pearson", use = "pairwise.complete") ** 2, 3)    
    print(r2)
    points(x = Ls, y = ybottom + scale * diff, col = col, type = "l")
    text(x = Ls[1], y = ybottom + scale * 1, labels = paste0("r2 = ", r2), pos = 1, cex = 1.5)
    return(r2)
}
add_numbers <- function(ytop, ybottom, x, i) {
    ## accept only if on a grid of every 50 SNPs
    ok <- unique(round(seq(1, length(x), length.out = 50)))
    ##q <- diff(range(L_grid)) / 20 ## bins
    w <- which((ytop - ybottom) > 0.001)
    w <- w[w %in% ok]
    if (length(w) > 0) {
        ww <- w
        ##a <- floor(x[w] / q)
        ##ww <- w[match(unique(a), a)]
        ## thin, make sure in different bins
        ## ww <- w[unique(round(seq(1, length(w), length.out = 20)))]
        text(x = x[ww], y = ((ytop + ybottom)/2)[ww], i)
    }
}
smooth_vector <- function(x, n) {
    x <- as.numeric(x)
    y <- (cumsum(x)[(n):(length(x))])
    z <- cumsum(c(0, x[1:(length(x) - n)]))
    a <- (y - z) / n
    return(a)
}
plot_single_gamma_dosage <- function(
    sampleReads,
    fbsoL,
    L_grid,
    L,
    inRegion2,
    cM_grid,
    ancAlleleFreqAll,
    haps,
    outname,
    method,
    truth_labels,
    have_truth_haplotypes,
    uncertain_truth_labels,
    output_plot = TRUE
) {
    ##
    colStore <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    smoothV <- 10
    Ls <- smooth_vector(L, smoothV)
    ##Ls <- smooth_vector(L_grid, smoothV)
    scale_dosage <- 0.5
    ##    
    nCols <- length(colStore)
    if (output_plot) {
        png(outname, height = 10, width = 20, units = "in", res = 100)
        par(mfrow = c(3, 1))
    }
    par(oma = c(0, 0, 5, 0))
    print(outname)
    ##
    ##
    ##
    r2s <- NULL
    for(i_which in 1:2) {
        par(mar = c(0, 0, 3, 0))
        if (method == "gibbs-nipt") {
            if (i_which == 1) { gammaK_t <- fbsoL$gammaMT_t }
            if (i_which == 2) { gammaK_t <- fbsoL$gammaMU_t }
            if (i_which == 3) { gammaK_t <- fbsoL$gammaP_t }
        } else if (method == "triploid-nipt") {
            if (i_which == 1) { gammaK_t <- fbsoL$list_of_gammas[[1]]$gammaKMT_t }
            if (i_which == 2) { gammaK_t <- fbsoL$list_of_gammas[[1]]$gammaKMU_t }
            if (i_which == 3) { gammaK_t <- fbsoL$list_of_gammas[[1]]$gammaKP_t }
        }
        ##
        K <- nrow(gammaK_t)
        xlim <- range(L_grid)
        ylim <- c(0, 1 + scale_dosage + scale_dosage)
        nGrids <- ncol(gammaK_t)
        backwards <- nGrids:1
        ##
        main <- c("hap1", "hap2")[i_which]
        ##
        plot(x = L_grid[1], y = 0, xlim = xlim, ylim = ylim, axes = FALSE, main = main, cex = 1.5)
        x <- L_grid ## c(L_grid[1], L_grid) ## , L_grid[length(L_grid):1])
        xleft <- c(x[1] - (x[2] - x[1]) / 2, x[-length(x)])
        xright <- c(x[-1], x[length(x)] + (x[length(x)] - x[(length(x) - 1)]) / 2)
        m <- array(0, c(nGrids, K + 1))
        ## is this slow...
        for(i in 1:K) {
            m[, i + 1] <- m[, i] + gammaK_t[i, ]
        }
        ##
        for(j in 1:2) {
            for(i in K:1) {
                ## can I o
                ybottom <- m[, i]
                ytop <- m[, i + 1]
                if (max(ytop - ybottom) > 0.01) {
                    if (j == 1) {
                        rect(
                            xleft = xleft,
                            xright = xright,
                            ybottom = ybottom,
                            ytop = ytop,
                            border = NA,
                            col = colStore[(i %% nCols) + 1]
                        )
                    } else {
                        add_numbers(ytop, ybottom, x, i)
                    }
                }
            }
        }
        ## distance to hap 1 and distance to hap2, same time?
        dosage <- fbsoL[["hapProbs_t"]][i_which, ]
        ## dosage2 <- fbsoL[["hapProbs_t"]][
        truth1 <- haps[, 1]
        truth2 <- haps[, 2]
        ## ## 
        r2s <- c(r2s, plot_1_dosage_vs_truth(dosage = dosage, truth = truth1, ancAlleleFreq = ancAlleleFreqAll, inRegion2 = inRegion2, smoothV = smoothV, Ls = Ls, ybottom = 1, scale = scale_dosage, col = "red"))
        r2s <- c(r2s, plot_1_dosage_vs_truth(dosage = dosage, truth = truth2, ancAlleleFreq = ancAlleleFreqAll, inRegion2 = inRegion2, smoothV = smoothV, Ls = Ls, ybottom = 1 + scale_dosage, scale = scale_dosage, col = "green"))
        ## add cor!
    }
    ## below, plot dosage vs truth
    recombF <- 0.5
    ylim <- c(0, 1 + recombF)
    plot(x = L_grid[1], y = 0, xlim = xlim, ylim = ylim, axes = FALSE, main = main, cex = 1.5)
    ## cM / Mbp difference
    rate <- diff(cM_grid) / diff(L_grid) * 1e6 ## normally 0-100 scaled
    rate <- rate / 100 * recombF
    rate[rate > recombF] <- recombF
    rate <- rate + 1
    lines(x = (L_grid[-1] + L_grid[-length(L_grid)]) / 2, y = rate, col = "purple", lwd = 2)
    abline(h = 1)
    ## add in how well the reads are doing
    test <- fbsoL$double_list_of_ending_read_labels[[1]][[1]]
    if (have_truth_haplotypes) {
        truth <- truth_labels
        truth[uncertain_truth_labels] <- 0
    }
    y <- 0.1 + (test - 1) / 2 + runif(length(test)) / 4
    ##
    for(iRead in 1:length(sampleReads)) {
        u <- range(sampleReads[[iRead]][[4]])
        ## level depending on what I say
        ## colour dpeending on truth
        if (have_truth_haplotypes) {
            col <- c("black", "blue", "red")[truth[iRead] + 1]
        } else {
            col <- "black"
        }
        lwd <- 1
        ## if (unhappy_reads[iRead]) {
        ##     col <- "purple"
        ##     lwd <- 5
        ## }
        segments(x0 = L[u[1] + 1], x1 = L[u[2] + 1], y0 = y[iRead], y1 = y[iRead], col = col, lwd = lwd)
    }
    ## 
    r2s <- c(r2s, plot_2_dosage_vs_truth(dosage = colSums(fbsoL[["hapProbs_t"]][1:2, ]), truth = rowSums(haps[, 1:2]), ancAlleleFreq = ancAlleleFreqAll, inRegion2 = inRegion2, smoothV = smoothV, Ls = Ls, ybottom = 0, scale = 1, col = "green"))
    outer_main <- "placeholder"
    mtext(text = outer_main, outer = TRUE, cex = 1.5)
    if (output_plot) {
        dev.off()
    }
    ## also add recombination rate map
    return(r2s)
}


robbie_image <- function(x, truth, n_breaks = 101, xlab = "", ylab = "") {
    ## 
    cols <- colorRampPalette(c("white", "red"))(n_breaks)
    x2A <- matrix(cols[cut(x, breaks = n_breaks)], nrow = nrow(x))
    cols <- colorRampPalette(c("white", "green"))(n_breaks)    
    x2B <- matrix(cols[cut(x, breaks = n_breaks)], nrow = nrow(x))
    ## 
    xlim <- c(0, ncol(x2A) + 1)
    ylim <- c(0, nrow(x2A) + 1)
    plot(x = 0, y= 0, xlim = xlim, ylim = ylim, col = "white", xlab = xlab, ylab = ylab)
    xleft <- 1:ncol(x2A) - 0.5
    xright <- xleft + 1
    ## 
    for(i_row in 1:nrow(x2A)) {
        col <- rep("white", ncol(x2A))
        w1 <- truth == 1 & !is.na(truth)
        col[w1] <- x2A[i_row, w1]
        w2 <- truth == 2 & !is.na(truth)
        col[w2] <- x2B[i_row, w2]
        ## run two gamuts, depending on match
        rect(xleft = xleft, xright = xright, ybottom = i_row - 0.5, ytop = i_row + 0.5, col = col, border = NA)
    }
}
