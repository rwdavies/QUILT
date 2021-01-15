plot_1_dosage_vs_truth <- function(dosage, truth, ancAlleleFreqAll, inRegion2, smoothV, Ls, scale, ybottom, col, minAFForBounding = 0.01, label = "") {
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
    points(x = Ls, y = ybottom + (0.8 * scale) * diff, col = col, type = "l")
    text_labels <- paste0("r2 = ", r2)
    if (label != "") {
        text_labels <- paste0(label, ", ", text_labels)
    }
    text(x = Ls[1], y = ybottom + scale * 0.85, labels = text_labels, pos = 4, cex = 1.25)
    return(r2)
}
plot_2_dosage_vs_truth <- function(dosage, truth, ancAlleleFreqAll, inRegion2, smoothV, Ls, scale, ybottom, col, minAFForBounding = 0.01, label = "") {
    ##
    ancAlleleFreqAll[ancAlleleFreqAll < minAFForBounding] <- minAFForBounding
    ancAlleleFreqAll[ancAlleleFreqAll > (1 - minAFForBounding)] <- (1 - minAFForBounding)
    ## 
    dosage2 <- dosage - 2 * ancAlleleFreqAll
    truth[is.na(truth)] <- ancAlleleFreqAll[is.na(truth)]
    truth2 <- truth - 2 * ancAlleleFreqAll
    ## difference!
    diff <- smooth_vector(x = abs(dosage2 - truth2), n = smoothV)
    ## diff <- diff # / 3
    diff[diff > 1] <- 1
    r2 <- round(cor(dosage2[inRegion2], truth2[inRegion2], method = "pearson", use = "pairwise.complete") ** 2, 3)
    points(x = Ls, y = ybottom + (0.8 * scale) * diff, col = col, type = "l")
    text_labels <- paste0("r2 = ", r2)
    if (label != "") {
        text_labels <- paste0(label, ", ", text_labels)
    }
    text(x = Ls[1], y = ybottom + scale * 0.85, labels = text_labels, pos = 4, cex = 1.25)
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
    sample_name,
    smooth_cm,
    regionStart,
    regionEnd,
    buffer,
    output_plot = TRUE
) {
    ##load("~/Downloads/atta/haps.NA12878HT.chr20.2000001.4000000_igs.1.it1.gibbs.png.RData")
   ## setwd("~/Downloads/atta")
    ##regionStart <- 1800000
    ##regionEnd <- 2000000
    ##buffer <- 200000
    ##outname="temp.png"
    ##
    colStore <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    smoothV <- 10
    Ls <- smooth_vector(L, smoothV)
    ##Ls <- smooth_vector(L_grid, smoothV)
    scale_dosage <- 0.5
    ##    
    nCols <- length(colStore)
    if (output_plot) {
        png(outname, height = 12, width = 20, units = "in", res = 100)
        par(mfrow = c(4, 1))
    }
    par(oma = c(0, 0, 5, 0))
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
        main <- c("Imputed haplotype 1", "Imputed haplotype 2")[i_which]
        ##
        plot(x = L_grid[1], y = 0, xlim = xlim, ylim = ylim, axes = FALSE, main = main, cex.main = 1.5, col = "white")
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
        r2s <- c(r2s, plot_1_dosage_vs_truth(dosage = dosage, truth = truth2, ancAlleleFreq = ancAlleleFreqAll, inRegion2 = inRegion2, smoothV = smoothV, Ls = Ls, ybottom = 1, scale = scale_dosage, col = "red", label = "Rolling accuracy versus truth haplotype 2"))
        r2s <- c(r2s, plot_1_dosage_vs_truth(dosage = dosage, truth = truth1, ancAlleleFreq = ancAlleleFreqAll, inRegion2 = inRegion2, smoothV = smoothV, Ls = Ls, ybottom = 1 + scale_dosage, scale = scale_dosage, col = "green", label = "Rolling accuracy versus truth haplotype 1"))
        ## add cor!
    }
    ##
    ## plot section 3 - add in how well the reads are doing, and the genotypes
    ##
    ylim <- c(0, 1)
    plot(x = L_grid[1], y = 0, xlim = xlim, ylim = ylim, axes = FALSE, main = "Read assignments", cex = 1.5, col = "white", cex.main = 1.5)
    test <- fbsoL$double_list_of_ending_read_labels[[1]][[1]] ## 1, 2, maybe 3 in the future
    if (have_truth_haplotypes) {
        truth <- truth_labels
        truth[uncertain_truth_labels] <- 0
    }
    ## so have space 0-1 to fit both of these
    ## use first two thirds for reads, last bit for diploid dosage
    frac_reads <- 2/3
    ## note, do 3-test so first haplotype reads on top
    y <- (1 - frac_reads) + (0.95 * frac_reads) * (((3 - test) - 1) * (2 / 3) + runif(length(test)) / 3)
    ##
    ## text(x = mean(xlim), y = (1 - frac_reads) + frac_reads / 2, labels = "Read assignments", cex = 1.5, font = 2)
    text(x = xlim[1], y = (1 - frac_reads) + frac_reads * (5 / 6), labels = "Hap1", srt = 90, cex = 1.5, font = 2)
    text(x = xlim[1], y = (1 - frac_reads) + frac_reads * (1 / 6), labels = "Hap2", srt = 90, cex = 1.5, font = 2)
    for(iRead in 1:length(sampleReads)) {
        u <- range(sampleReads[[iRead]][[4]])
        ## level depending on what I say
        ## colour dpeending on truth
        if (have_truth_haplotypes) {
            col <- c("black", "blue", "red")[truth[iRead] + 1]
            text(x = Ls[1], y = (1 - frac_reads) + frac_reads / 2, labels = "Blue = truth hap 1, Red = truth hap 2", pos = 4, cex = 1.25)
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
    ## abline(h = (1 - frac_reads))
    text(x = mean(xlim), y = (1 - frac_reads) * 0.85, labels = "Imputed genotypes", cex = 1.5, font = 2)
    r2s <- c(r2s, plot_2_dosage_vs_truth(dosage = colSums(fbsoL[["hapProbs_t"]][1:2, ]), truth = rowSums(haps[, 1:2]), ancAlleleFreq = ancAlleleFreqAll, inRegion2 = inRegion2, smoothV = smoothV, Ls = Ls, ybottom = 0, scale = (1 - frac_reads), col = "green", label = "Rolling accuracy versus truth genotypes"))
    ##
    ## plot section 4 - recombination rate (need real thing!), and x-axis with locations
    ##
    recombF <- 0.3 ## i.e. what fraction of the below is for recombination
    ylim <- c(0, recombF)
    ## the extra 1 space is for the rest below
    ##rate <- diff(cM_grid) / diff(L_grid) * 1e6 ## normally 0-100 scaled
    rate <- smooth_cm / max(smooth_cm) ## 0-1 scale
    rate <- rate * recombF ## so on a scale from 0 to recombF
    ## recombF <- max(rate) ## bound slightly?
    par(mar = c(5, 0, 3, 0)) ## bottom, left, top, and right.
    plot(x = L_grid[1], y = 0, xlim = xlim, ylim = ylim, axes = FALSE, main = "Smoothed recombination rate", cex = 1.5, col = "white", cex.main = 1.5, xlab = "Physical position")
    lines(x = (L_grid[-1] + L_grid[-length(L_grid)]) / 2, y = rate, col = "purple", lwd = 2)
    axis(1)
    ##
    ## finalize
    ##
    outer_main <- paste0("QUILT imputation for ", sample_name)
    mtext(text = outer_main, outer = TRUE, cex = 1.5)
    abline(v = regionStart)
    if (!is.na(regionStart)) {
        text(x = regionStart, srt = 90, labels = "Region Start", cex = 1.5, font = 2, pos = 1)
        text(x = regionEnd, srt = 90, labels = "Region End", cex = 1.5, font = 2, pos = 3)
    }
    abline(v = regionStart)
    abline(v = regionEnd)
    if (output_plot) {
        dev.off()
    }
    ## abline(h = 1)
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
