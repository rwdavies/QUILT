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
add_numbers <- function(ytop, ybottom, x, i, col = "black", n = 50) {
    ## accept only if on a grid of every 50 SNPs
    ok <- unique(round(seq(1, length(x), length.out = n)))
    ##q <- diff(range(L_grid)) / 20 ## bins
    w <- which((ytop - ybottom) > 0.001)
    w <- w[w %in% ok]
    if (length(w) > 0) {
        ww <- w
        ##a <- floor(x[w] / q)
        ##ww <- w[match(unique(a), a)]
        ## thin, make sure in different bins
        ## ww <- w[unique(round(seq(1, length(w), length.out = 20)))]
        text(x = x[ww], y = ((ytop + ybottom)/2)[ww], i, col = col)
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
    have_truth_genotypes,
    truth_gen,
    uncertain_truth_labels,
    sample_name,
    smooth_cm,
    regionStart,
    regionEnd,
    buffer,
    new_haps = NULL,
    which_haps_to_use = NULL,
    output_plot = TRUE
) {
    ##
    colStore <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    smoothV <- 10
    Ls <- smooth_vector(L, smoothV)
    ##Ls <- smooth_vector(L_grid, smoothV)
    scale_dosage <- 0.5
    n <- c(diploid = 2, nipt = 3)[method]
    stopifnot(!is.na(n))
    ##
    nCols <- length(colStore)
    if (output_plot) {
        png(outname, height = (n + 3) * 3, width = 20, units = "in", res = 100)
        par(mfrow = c(n + 3, 1))
    }
    par(oma = c(0, 0, 5, 0))
    ##
    ##
    ##
    r2s <- NULL
    for(i_which in 1:n) {
        par(mar = c(0, 0, 3, 0))
        if (i_which == 1) { gammaK_t <- fbsoL$gammaMT_t }
        if (i_which == 2) { gammaK_t <- fbsoL$gammaMU_t }
        if (i_which == 3) { gammaK_t <- fbsoL$gammaP_t }
        ## } else if (method == "triploid-nipt") {
        ##     if (i_which == 1) { gammaK_t <- fbsoL$list_of_gammas[[1]]$gammaKMT_t }
        ##     if (i_which == 2) { gammaK_t <- fbsoL$list_of_gammas[[1]]$gammaKMU_t }
        ##     if (i_which == 3) { gammaK_t <- fbsoL$list_of_gammas[[1]]$gammaKP_t }
        ## }
        ##
        K <- nrow(gammaK_t)
        if (length(new_haps) > 0) {
            letter_col <- rep("black", K)
            letter_col[new_haps] <- "white"
        } else {
            letter_col <- rep("black", K)
        }
        xlim <- range(L_grid)
        ylim <- c(0, 1 + scale_dosage + scale_dosage)
        if (method == "nipt") {
            ylim[2] <- ylim[2] + scale_dosage
        }
        nGrids <- ncol(gammaK_t)
        backwards <- nGrids:1
        ##
        main <- paste0("Imputed haplotype ", i_which)
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
        if (is.null(which_haps_to_use)) {
            which_haps_to_use <- 1:K
        }
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
                        ## add_numbers(ytop, ybottom, x, i, col = letter_col[i])
                        add_numbers(ytop, ybottom, x, which_haps_to_use[i]) ## make global
                    }
                }
            }
        }
        ## distance to hap 1 and distance to hap2, same time?
        dosage <- fbsoL[["hapProbs_t"]][i_which, ]
        ## dosage2 <- fbsoL[["hapProbs_t"]][
        if (have_truth_haplotypes) {
            if (method == "nipt") {
                truth3 <- haps[, 3]
                r2s <- c(r2s, plot_1_dosage_vs_truth(dosage = dosage, truth = truth3, ancAlleleFreq = ancAlleleFreqAll, inRegion2 = inRegion2, smoothV = smoothV, Ls = Ls, ybottom = 1, scale = scale_dosage, col = "red", label = "Rolling accuracy versus truth haplotype 3"))
                offset <- scale_dosage
            } else {
                offset <- 0
            }
            truth1 <- haps[, 1]
            truth2 <- haps[, 2]
            r2s <- c(r2s, plot_1_dosage_vs_truth(dosage = dosage, truth = truth2, ancAlleleFreq = ancAlleleFreqAll, inRegion2 = inRegion2, smoothV = smoothV, Ls = Ls, ybottom = 1 + offset, scale = scale_dosage, col = "red", label = "Rolling accuracy versus truth haplotype 2"))
            r2s <- c(r2s, plot_1_dosage_vs_truth(dosage = dosage, truth = truth1, ancAlleleFreq = ancAlleleFreqAll, inRegion2 = inRegion2, smoothV = smoothV, Ls = Ls, ybottom = 1 + offset + scale_dosage, scale = scale_dosage, col = "green", label = "Rolling accuracy versus truth haplotype 1"))
        }
        ## add cor!
    }
    ##
    ## plot section 3 - add in how well the reads are doing, and the genotypes
    ##
    ylim <- c(0, 1)
    plot(x = L_grid[1], y = 0, xlim = xlim, ylim = ylim, axes = FALSE, main = "Read assignments", cex = 1.5, col = "white", cex.main = 1.5, xlab= "", ylab = "")
    test <- fbsoL$double_list_of_ending_read_labels[[1]][[1]] ## 1, 2, maybe 3 in the future
    if (have_truth_haplotypes) {
        truth <- truth_labels
        truth[uncertain_truth_labels] <- 0
    }
    ## so have space 0-1 to fit both of these (or three)
    ## use first two thirds for reads, last bit for diploid dosage
    if (method == "diploid") {
        frac_reads <- 2/3
    } else{
        frac_reads <- 3 / 4
    }
    ## note, do 3-test so first haplotype reads on top
    ##
    ## text(x = mean(xlim), y = (1 - frac_reads) + frac_reads / 2, labels = "Read assignments", cex = 1.5, font = 2)
    if (method == "diploid") {
        y <- (1 - frac_reads) + (0.95 * frac_reads) * (((3 - test) - 1) * (2 / 3) + runif(length(test)) / 3)
        text(x = xlim[1], y = (1 - frac_reads) + frac_reads * (5 / 6), labels = "Hap1", srt = 90, cex = 1.5, font = 2)
        text(x = xlim[1], y = (1 - frac_reads) + frac_reads * (1 / 6), labels = "Hap2", srt = 90, cex = 1.5, font = 2)
    } else {
        ## put these at the top
        y <- rep((1 - frac_reads), length(test))
        H <- test
        y[H == 1] <- y[H == 1] + (frac_reads) * (4 / 5 + runif(sum(H == 1)) * 1 / 5)
        y[H == 2] <- y[H == 2] + (frac_reads) * (2 / 5 + runif(sum(H == 2)) * 1 / 5)
        y[H == 3] <- y[H == 3] + (frac_reads) * (0 / 5 + runif(sum(H == 3)) * 1 / 5)
        text(x = xlim[1], y = (1 - frac_reads) + frac_reads * (9 / 10), labels = "Hap1", srt = 90, cex = 1.5, font = 2)
        text(x = xlim[1], y = (1 - frac_reads) + frac_reads * (5 / 10), labels = "Hap2", srt = 90, cex = 1.5, font = 2)
        text(x = xlim[1], y = (1 - frac_reads) + frac_reads * (1 / 10), labels = "Hap3", srt = 90, cex = 1.5, font = 2)
    }
    ##
    if (have_truth_haplotypes) {
        if (method == "nipt") {
            label <- "Orange = truth hap 1, Green = truth hap 2, Purple = truth hap 3"
        } else {
            label <- "Orange = truth hap 1, Green = truth hap 2"
        }
        text(
            x = Ls[1], y = (1 - frac_reads) + frac_reads / 2,
            labels = label,
            pos = 4, cex = 1.25
        )
    }
    us <- t(sapply(1:length(sampleReads), function(iRead) {
        range(sampleReads[[iRead]][[4]])
    }))
    ## level depending on what I say
    ## colour dpeending on truth
    grey <- alpha_col("grey", 0.25)
    if (have_truth_haplotypes) {
        col <- c(grey, "orange", "green", "purple")[truth + 1]
    } else {
        col <- rep(grey, length(test))
    }
    lwd <- 1
    ## over-plot coloured ones
    for(i in 1:2) {
        if (i == 1) w <- col == grey
        if (i == 2) w <- col != grey
        segments(x0 = L[us[w, 1] + 1], x1 = L[us[w, 2] + 1], y0 = y[w], y1 = y[w], col = col[w], lwd = lwd)
    }
    ##
    ## abline(h = (1 - frac_reads))
    text(x = mean(xlim), y = (1 - frac_reads) * 0.85, labels = "Imputed genotypes", cex = 1.5, font = 2)
    ##
    ##
    ##
    if (have_truth_haplotypes | have_truth_genotypes) {
        if (method == "nipt") {
            if (have_truth_haplotypes) {
                truthH_m <- rowSums(haps[, 1:2])
                truthH_g <- rowSums(haps[, c(1, 3)])
                ## not sure what to do about truth_genotypes here honestly
                r2s <- c(r2s, plot_2_dosage_vs_truth(dosage = colSums(fbsoL[["hapProbs_t"]][1:2, ]), truth = truthH_m, ancAlleleFreq = ancAlleleFreqAll, inRegion2 = inRegion2, smoothV = smoothV, Ls = Ls, ybottom = (1 - frac_reads) / 2, scale = (1 - frac_reads) / 2, col = "green", label = "Rolling accuracy of maternal genotypes versus truth"))
                r2s <- c(r2s, plot_2_dosage_vs_truth(dosage = colSums(fbsoL[["hapProbs_t"]][c(1, 3), ]), truth = truthH_g, ancAlleleFreq = ancAlleleFreqAll, inRegion2 = inRegion2, smoothV = smoothV, Ls = Ls, ybottom = 0, scale = (1 - frac_reads) / 2, col = "orange", label = "Rolling accuracy of fetal genotypes versus truth"))
            }
        } else {
            if (have_truth_genotypes & !have_truth_haplotypes) {
                truthH <- c(truth_gen)
            } else {
                truthH <- rowSums(haps[, 1:2])
            }
            r2s <- c(r2s, plot_2_dosage_vs_truth(dosage = colSums(fbsoL[["hapProbs_t"]][1:2, ]), truth = truthH, ancAlleleFreq = ancAlleleFreqAll, inRegion2 = inRegion2, smoothV = smoothV, Ls = Ls, ybottom = 0, scale = (1 - frac_reads), col = "green", label = "Rolling accuracy versus truth genotypes"))
        }
    }
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





## Modified from "benchnmark_functions"
plot_of_likelihood_with_time_new <- function(
    p_store,
    n_hgfi = 1,
    n_samplings = 1,
    ylim_options = NULL,
    option = NULL,
    outname = NULL,
    y_to_plot = c("p_O_given_H_L", "p_H_given_L", "p_H_given_O_L_up_to_C"),
    truth_likelihoods = NULL
) {
    if (!is.null(option)) {
        outer_main <- paste0(names(option), as.character(option), sep = "=", collapse = ", ")
    } else {
        outer_main <- NULL
    }
    if (is.null(ylim_options)) {
        ylim_options <- c("bounded")
    }
    ## hmm, what to plot...

    for(ylim_option in ylim_options) {
        if (!is.null(outname)) {
            png(outname, height = length(y_to_plot) * 5, width = 20, units = "in", res = 300)
        }
        par(mfrow = c(3, 1))
        cbPalette <- rep(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), ceiling(n_samplings / 8))
        ##
        for(ycol in y_to_plot) {
            if (ylim_option == "bounded") {
                ylim <- c(-100, 0)
                ylab <- "log10-likelihood + C"
            } else if (ylim_option == "bounded5000") {
                ylim <- c(-5000, 0)
                ylab <- "log10-likelihood + C"
            } else if (ylim_option == "full") {
                ylim <- range(p_store[, ycol])
                ylab <- "log10-likelihood"
            } else {
                stop("bad ylim_option")
            }
            ##
            nReads <- nrow(p_store) / n_samplings / n_hgfi
            x <- 1:(nReads * n_hgfi)
            xlim <- c(1, nReads * n_hgfi)
            plot(
                x = 0, y = 0, xlim = xlim, ylim = ylim,
                xlab = "Increasing parameter updates ->\nRed line = Complete update of all reads",
                ylab = ylab,
                axes = FALSE,
                col = "white",
                main = ycol
            )
            for(i in 0:(n_hgfi + 1)) {
                abline(v = i * nReads, col = "red")
            }
            ## p1, p2, p3, p, p_H_given_L
            ## what about pSR? etc?
            if (ylim_option == "bounded") {
                t <- p_store[, ycol]
                ## if more than 2 iterations, get max from later!
                if (n_hgfi > 2) {
                    ## if initializing iteratively, take the max after the first two iterations, if possible
                    xx <- 1:(nReads * 2)
                    for(i_sampling in 1:n_samplings) {
                        t[(i_sampling - 1) * n_hgfi * nReads + xx] <- NA
                    }
                }
                max_val <- max(t, na.rm = TRUE)
            }
            ## add truth here
            if (!is.null(truth_likelihoods)) {
                for(i_ll in 1:nrow(truth_likelihoods)) {
                    prob <- truth_likelihoods[i_ll, ycol]
                    if (ylim_option == "bounded") {
                        prob <- prob - max_val
                    }
                    abline(h = prob, col = "green")
                }
            }
            ##
            for(i_sampling in 1:n_samplings) {
                probs <- p_store[
                (i_sampling - 1) * n_hgfi * nReads + x, ycol
                ]
                if (ylim_option == "bounded") {
                    probs <- probs - max_val
                }
                points(probs, col = cbPalette[i_sampling], type = "l")
            }

            ##axis(1)
            axis(2)
        }
        if (!is.null(outname)) {
            dev.off()
        }
    }
}





## Modified from "benchnmark_functions"
plot_of_likelihoods_across_samplings_and_seek_its <- function(
    outputdir,
    for_likelihood_plotting,
    nGibbsSamples,
    n_gibbs_burn_in_its,
    block_gibbs_iterations,
    n_seek_its,
    sample_name,
    regionName,
    truth_likelihoods = NULL,
    ylab = "p_H_given_O_L_up_to_C"
) {
    cbPalette <- rep(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), ceiling(nGibbsSamples / 8))
    ## each of nGibbssSamples is a different colour
    ## n_seek_its is blocks
    ##
    ## initialize
    ##
    filename <- file.path(
        outputdir,
        "plots",
        paste0("likelihoods.", sample_name, ".",regionName, ".png")
    )
    png(
        filename,
        height = 2000,
        width = 2000 * n_seek_its * (n_gibbs_burn_in_its / 10),
        res = 300
    )
    par(mfrow = c(2, 1))
    for(ylab_option in c("zoomedout", "zoomedin")) {
        ## first option is zoomed out
        ## second option is zoomed in
        ##
        ## determine what ylimit should be?
        ##
        xlim <- c(0, 1)
        y <- sapply(for_likelihood_plotting, function(x) {
            sapply(x, function(y) {
                if (ylab_option == "zoomedout") {
                    if (n_gibbs_burn_in_its > 4) {
                        y[-c(1:2), ylab]
                    }
                    y[, ylab]
                } else {
                    ## just keep last 20% (or fewer)
                    y[max(1, round(0.20 * nrow(y)) - 3):nrow(y), ylab]
                }
            })
        })
        ylim <- range(y)
        plot(
            x = 0, y = 0,
            xlim = xlim, ylim = ylim,
            col = "white",
            xlab = "Progression through Gibbs sampling",
            ylab = "- log likelihood + constant",
            axes = FALSE
        )
        axis(2)
        ## plot seek its
        abline(v = 0, col = "black")
        for(i_seek_it in 1:n_seek_its) {
            abline(v = i_seek_it / n_seek_its, col = "black")
            ## plot block Gibbs
            x <- seq(
            (i_seek_it - 1) / n_seek_its,
            (i_seek_it - 0)/ n_seek_its,
            length.out = nrow(for_likelihood_plotting[[1]][[1]])
            )
            if (length(block_gibbs_iterations) > 0) {
                for(val in x[block_gibbs_iterations]) {
                    abline(v = val, col = "purple")
                }
            }
        }
        for(iGibbsSample in 1:nGibbsSamples) {
            for(i_seek_it in 1:n_seek_its) {
                data <- for_likelihood_plotting[[iGibbsSample]][[i_seek_it]]
                y <- data[, ylab]
                col <- cbPalette[iGibbsSample]
                x <- seq(
                (i_seek_it - 1) / n_seek_its,
                (i_seek_it - 0)/ n_seek_its,
                length.out = length(y)
                )
                points(x = x, y = y, type = "l", col = col)
            }
            ##
        }
    }
    dev.off()
    ## system(paste("rsync -av ", filename, " rescompNew2:~/"))
}




plot_prob_of_flipping_to_first_hap <- function(
    p_store,
    filename
) {

    filename <- "~/temp.png"
    load(file = "~/temp.RData")

    nReads <- ncol(p1_store[[1]][[1]])
    n_sampling_its <- nrow(p1_store[[1]][[1]])

    n_seek_its <- length(p1_store[[1]])
    nGibbsSamples <- length(p1_store) - 1

    ##
    size <- round(nReads / 1000)
    size <- min(max(size, 5), 20)

    ##
    ## here look at how probabilities of sampling change
    ##
    for(i_what in 1:2) {
        if (i_what == 1) {
            filename <- "~/temp.readprobs.png"
        } else {
            filename <- "~/temp.reads.png"
        }
        png(filename, height = size, width = size, units = "in", res = 300)
        xlim <- c(1, n_seek_its + 1)
        ylim <- c(1, nGibbsSamples + 1 + 1)
        plot(x = 0, y = 0, col = "white", xlim = xlim, ylim = ylim, axes = TRUE, xlab = "", ylab = "")
        ybottom <- (0:(nReads - 1)) / nReads
        ytop <- ybottom + 1 / nReads
        xleft <- 0:(n_sampling_its - 1) / n_sampling_its
        xright <- xleft + 1 / n_sampling_its
        ## now rep them
        ybottom <- rep(ybottom, each = n_sampling_its)
        ytop <- rep(ytop, each = n_sampling_its)
        ##
        xleft <- rep(xleft, nReads)
        xright <- rep(xright, nReads)
        ##
        cols <- colorRampPalette(c("#FFFFFF", "#56B4E9"))(51)
        for(i_it in 1:n_seek_its) {
            for(iGibbs in 1:(nGibbsSamples + 1)) {
                if (i_what == 1) {  p1 <- p1_store[[iGibbs]][[i_it]] }
                if (i_what == 2) {  p1 <- read_store[[iGibbs]][[i_it]] }
                rect(i_it + xleft, iGibbs + ybottom, i_it + xright, iGibbs + ytop, col = cols[round(100 * abs(p1 - 0.5)) + 1], border = NA)
                rect(i_it, iGibbs, i_it + 1, iGibbs + 1)
            }
        }
        dev.off()
    }


    ##
    ## here plot the reads themselves
    ##
    table(round(colMeans(abs(x - 0.5)), 1))

    ## plot reads used as well
    filename <- "temp.reads.png"
    png(filename, height = size, width = size, units = "in", res = 300)

    ## argh, not yet recorded
    f <- function(i, j) {
        x <- p1_store[[i]][[j]]
        y <- round(apply(x, 2, function(x) max(abs(x - 0.5))), 1)
        y
    }

    table(f(1, 1), f(1, 3))
    table(f(2, 3), f(1, 3)) ## mostly the same


    ##
    ## here compare the read probs at the start and end of each run
    ##
    o <- lapply(p1_store, function(p1) {
        m <- nrow(p1_store[[1]][[1]])
        m <- cbind(
            t(p1[[1]][c(1, m), ]),
            t(p1[[2]][c(1, m), ]),
            t(p1[[3]][c(1, m), ])
        )
        m
    })
    a <- cbind(o[[1]], o[[2]])
    for(j in 3:length(o)) {
        a <- cbind(a, o[[j]])
    }
    a <- t(a)

    pdf("~/temp.pdf", height = 10, width = 30)
    xlim <- c(0, 1)
    ylim <- c(0, 1)
    plot(x = 0, y = 0, col = "white", xlim = xlim, ylim = ylim, axes = TRUE, xlab = "", ylab = "")
    ##
    ybottom <- (0:(nReads - 1)) / nReads
    ytop <- ybottom + 1 / nReads
    n_sampling_its <- nrow(a)
    xleft <- 0:(n_sampling_its - 1) / n_sampling_its
    xright <- xleft + 1 / n_sampling_its
    ## now rep them
    ybottom <- rep(ybottom, each = n_sampling_its)
    ytop <- rep(ytop, each = n_sampling_its)
    xleft <- rep(xleft, nReads)
    xright <- rep(xright, nReads)
    ##
    rect(ybottom, xleft, ytop, xright, col = cols[round(100 * abs(a - 0.5)) + 1], border = NA)
    ## add some lines
    abline(h = seq(0, n_sampling_its - 1, 6) / n_sampling_its)
    dev.off()

}



## this can be used to plot H too
plot_H_class <- function(H_class, L_grid, wif0) {
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    xlim <- range(L_grid)
    x <- L_grid
    a <- (x[-1] + x[-length(x)]) / 2
    xleft <- c(x[1] - (x[2] - x[1]) / 2, a)
    xright <- c(a, x[length(x)] + (x[length(x)] - x[length(x) - 1]) / 2)
    f8 <- function(x) {
        sapply(1:6, function(i) {
            sum(x == i)
        })
    }
    m1 <- t(sapply(as.integer(names(table(wif0))), function(iGrid) {
        w <- wif0 == iGrid ## meh
        f8(H_class[w])
    }))
    ## skip 0 and 7
    m1b <- t(apply(m1, 1, function(x) x / sum(x)))
    plot(x = 0, y = 0, xlim = xlim, ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
    ybottom <-rep(0, nrow(m1b))
    ytop <- rep(0, nrow(m1b))
    n <- ncol(m1b)
    for(i in 1:n) {
        ytop <- ytop + m1b[, i]
        rect(xleft = xleft, xright = xright, ybottom = ybottom, ytop = ytop, col = cbPalette[i], border = NA)
        ybottom <- ytop
    }
    NULL
}

plot_log_H_class_in_grid <- function(H_class, L_grid, wif0, ff, plot = TRUE, plot_new = FALSE, col = "black") {
    xlim <- range(L_grid)
    x <- L_grid
    a <- (x[-1] + x[-length(x)]) / 2
    xleft <- c(x[1] - (x[2] - x[1]) / 2, a)
    xright <- c(a, x[length(x)] + (x[length(x)] - x[length(x) - 1]) / 2)
    a <- log(c(0.5, 0.5 - ff / 2, ff / 2))
    m1 <- t(sapply(as.integer(names(table(wif0))), function(iGrid) {
        w <- wif0 == iGrid ## meh
        get_log_p_H_class(H_class[w], ff)  / sum(w)
    }))
    ## skip 0 and 7
    ylim <- c(min(a), 0) ## meh, should be OK
    if (plot) {
        if (plot_new) {
            plot(x = 0, y = 0, xlim = xlim, ylim = ylim, xlab = "", ylab = "", axes = FALSE)
        }
        lines(x = (xleft + xright) / 2, y = m1, lwd = 2, col = col)
    }
    return(m1)
}









##
## mostly just for interactive debugging
##
plot_shard_block_output <- function(
    shard_block_results,
    before_gamma1_t,
    before_gamma2_t,
    before_gamma3_t,
    after_gamma1_t,
    after_gamma2_t,
    after_gamma3_t,
    before_read_labels,
    after_read_labels,
    nGrids,
    outname,
    L_grid,
    L,
    uncertain_truth_labels,
    truth_labels,
    have_truth_haplotypes,
    sampleReads,
    wif0,
    grid,
    method,
    ff,
    only_plot_confident_reads,
    before_H_class,
    after_H_class,
    is_block = FALSE
) {
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    cbPalette4 <- sapply(cbPalette, alpha_col, alpha = 0.4)
    if (method != "nipt") {
        stop("function is partially adapted, but only for NIPT")
    }
    ##
    ##
    ##
    s <- 1
    S <- 1
    ##
    xlim <- range(L_grid)
    ##
    ##
    ##
    ##
    ## width <- min(max(20, (L_grid[length(L_grid)] - L_grid[1]) / 1e6 * 36), 200)
    height <- 20
    png(outname, height = height, width = 25, res = 200, units = "in")
    mfrow <- c(diploid = 11, nipt = 15)[method]
    par(mfrow = c(mfrow, 1))
    par(oma = c(0, 0, 5, 0))
    grid_distances <- diff(L_grid)
    x <- L_grid[-1] - grid_distances
    xleft <- L_grid[-length(L_grid)]
    xright <- L_grid[-1]
    ## for fbdstore
    midpoints <- L_grid[-1] - grid_distances / 2
    xleft2 <- c(L_grid[1], midpoints)
    xright2 <- c(midpoints, L_grid[length(L_grid)])
    ##
    ## define (for gamma), xleft and xright as follows
    ## L_grid will be the middle of the two
    ## x_left and x_right will average to L_grid
    ## with a first x_right and the next x_left adding up
    ##
    x <- L_grid
    a <- (x[-1] + x[-length(x)]) / 2
    xleft <- c(x[1] - (x[2] - x[1]) / 2, a)
    xright <- c(a, x[length(x)] + (x[length(x)] - x[length(x) - 1]) / 2)
    uu <- sapply(sampleReads, function(x) range(x[[4]])) + 1
    nReads <- length(sampleReads)
    for(i_type in 1:2) {
        ##
        ##
        ##
        if (i_type == 1) {
            read_labels <- before_read_labels
            gamma1_t <- before_gamma1_t
            gamma2_t <- before_gamma2_t
            gamma3_t <- before_gamma3_t
            what_we_are_plotting <- "before"
            H_class <- before_H_class
        } else {
            read_labels <- after_read_labels
            gamma1_t <- after_gamma1_t
            gamma2_t <- after_gamma2_t
            gamma3_t <- after_gamma3_t
            what_we_are_plotting <- "after"
            H_class <- after_H_class            
        }
        ##
        ## if this is the first one, add the comparison of H_class results
        if (i_type == 1) {
            par(mar = c(0, 0, 3, 0))
            m1 <- plot_log_H_class_in_grid(before_H_class, L_grid, wif0, ff, plot = TRUE, plot_new = TRUE, col = "black")
            m2 <- plot_log_H_class_in_grid(after_H_class, L_grid, wif0, ff, plot = TRUE, plot_new = FALSE, col = "red")
        }
        ## m2 <- plot_log_H_class_in_grid(after_H_class, L_grid, wif0, ff, plot = )
        ## par(mar = c(0, 0, 1, 0))
        ## xlim <- range(L_grid)
        ## x <- L_grid
        ## a <- (x[-1] + x[-length(x)]) / 2
        ## xleft <- c(x[1] - (x[2] - x[1]) / 2, a)
        ## xright <- c(a, x[length(x)] + (x[length(x)] - x[length(x) - 1]) / 2)
        ## plot(x = (xleft + xright) / 2, y = as.numeric(m2 / m1), type = "l", axes = FALSE) ## lower = better
        ## abline(h = 1)
        ##
        ## par(mar = c(0, 0, 3, 0))
        ## ylim <- c(0, 1)
        ## plot(x = L_grid[1], y = 0, xlim = xlim, ylim = ylim, axes = FALSE, cex = 1.5, col = "white")
        ## y <- 0.1 + (2 - (read_labels - 1)) / 3 + runif(length(read_labels)) / 4
        ## col <- cbPalette[before_H_class + 1]
        ## segments(x0 = L[uu[1, ]], x1 = L[uu[2, ]], y0 = y, y1 = y, col = col, lwd = 1.5)
        ## legend("topright", c("None", "1", "2", "3", "12", "13", "23", "123"), col = cbPalette, lwd = 2, cex = 0.75)
        par(mar = c(0, 0, 3, 0))
        if (i_type == 1) {
            H_class <- before_H_class
        } else {
            H_class <- after_H_class
        }
        plot_H_class(H_class, L_grid, wif0)
        ##
        ## if this is the second one, add the results in here
        ##
        if (i_type == 1) {
            ir_chosen <- rep(1, nrow(shard_block_results))
            for(jjj in 1:2) {
                ## first do left
                ## then do right
                if (jjj == 1) {
                    ir_col <- "ir_left"
                    prefix <- "left_p"
                } else {
                    ir_col <- "ir_right"
                    prefix <- "p"
                }
                ir_chosen_local <- shard_block_results[, ir_col]
                ir_chosen[ir_chosen_local != 1] <- 2 ## doesn't matter what number
                par(mar = c(0, 0, 3, 0))
                ylim <- c(0, 1.2)
                plot(x = L_grid[1], y = 0, xlim = xlim, ylim = ylim, axes = FALSE, cex = 1.5, col = "white")
                col <- cbPalette[ir_chosen_local]
                rect(xleft = xleft[-1], xright = xright[-1], ybottom = 1, ytop = 1.2, col = col, border = NA)
                ##
                ## ybottom <- rep(0, nrow(shard_block_results))
                ## ytop <- shard_block_results[, paste0(prefix, "1")]
                ## for(ir in 1:6) {
                ##     rect(xleft = xleft[-1], xright = xright[-1], ybottom = ybottom, ytop = ytop, col = cbPalette[ir], border = NA)
                ##     ybottom <- ytop
                ##     if (ir < 6) {
                ##         ytop <- ytop + shard_block_results[, paste0(prefix, ir + 1)]
                ##     }
                ## }
                abline(h = 1, lwd = 2)
                ## add vertical bars, for sure!
                changes <- which(diff(ir_chosen_local) != 0) + 1
                xa <- (xleft[-1])[changes]
                abline(v = xa, col = "black")
            }
            ## for the rest, add the probs
            ##changes <- which(diff(ir_chosen) != 0) + 1
            ##xa <- (xleft[-1])[changes]
            ##abline(v = xa, col = "black")
            ##if (length(changes) > 0) {
            ##    text(x = xa, y = 1, labels = ir_chosen[changes], col = "black")
            ##}
        }
        ##
        ## 2) add in reads here
        ##
        par(mar = c(0, 0, 3, 0))
        ylim <- c(0, 1)
        plot(x = L_grid[1], y = 0, xlim = xlim, ylim = ylim, axes = FALSE, cex = 1.5, col = "white")
        if (have_truth_haplotypes) {
            truth <- truth_labels
            truth[uncertain_truth_labels] <- 0
            if (method == "nipt") {
                label <- "Orange = truth hap 1, Blue = truth hap 2, Green = truth hap 3"
            } else {
                label <- "Orange = truth hap 1, Green = truth hap 2"
            }
            text(
                x = L[1], y = (1),
                labels = label,
                pos = 4, cex = 1.25, xpd = NA
            )
        }
        if (method == "nipt") {
            y <- 0.1 + (2 - (read_labels - 1)) / 3 + runif(length(read_labels)) / 4
        } else {
            y <- 0.1 + (read_labels - 1) / 3 + runif(length(read_labels)) / 4
        }
        lwd <- 1.5
        ##
        ##
        if (have_truth_haplotypes) {
            col <- cbPalette[truth + 1]
        } else {
            col <- rep(cbPalette[1], nReads)
        }
        if (only_plot_confident_reads & have_truth_haplotypes) {
            w <- truth != 0
            segments(x0 = L[uu[1, !w]], x1 = L[uu[2, !w]], y0 = y[!w], y1 = y[!w], col = alpha_col("grey", 0.4), lwd = lwd)
        } else {
            w <- rep(TRUE, nReads)
        }
        segments(x0 = L[uu[1, w]], x1 = L[uu[2, w]], y0 = y[w], y1 = y[w], col = col[w], lwd = lwd)
        ##
        ##
        ## 3) plot gammas
        ##
        ## for now, plot gammas before
        ##
        n <- c(diploid = 2, nipt = 3)[method]
        for(i_which in 1:n) {
            par(mar = c(0, 0, 3, 0))
            scale_dosage <- 0
            colStore <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
            nCols <- length(colStore)
            if (i_which == 1) { gammaK_t <- gamma1_t;   main <- paste0("Hap 1 - ", what_we_are_plotting)}
            if (i_which == 2) { gammaK_t <- gamma2_t;   main <- paste0("Hap 2 - ", what_we_are_plotting)}
            if (i_which == 3) { gammaK_t <- gamma3_t;   main <- paste0("Hap 3 - ", what_we_are_plotting)}
            ##
            K <- nrow(gammaK_t)
            ylim <- c(0, 1 + c(diploid = 2, nipt = 3)[method] * scale_dosage)
            nGrids <- ncol(gammaK_t)
            backwards <- nGrids:1
            ##
            plot(x = L_grid[1], y = 0, xlim = xlim, ylim = ylim, axes = FALSE, cex = 1.5, main = main)
            x <- L_grid ## c(L_grid[1], L_grid) ## , L_grid[length(L_grid):1])
            ##xleft <- c(x[1] - (x[2] - x[1]) / 2, x[-length(x)])
            ##xright <- c(x[-1], x[length(x)] + (x[length(x)] - x[(length(x) - 1)]) / 2)
            a <- (x[-1] + x[-length(x)]) / 2
            xleft <- c(x[1] - (x[2] - x[1]) / 2, a)
            xright <- c(a, x[length(x)] + (x[length(x)] - x[length(x) - 1]) / 2)
            ##
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
                            add_numbers(ytop, ybottom, x, i, n = round(nGrids / 5))
                        }
                    }
                }
            }
            ## always add in vertical bars at switch points
            changes <- which(diff(ir_chosen) != 0) + 1
            xa <- (xleft[-1])[changes]
            abline(v = xa, col = "black")
            if (length(changes) > 0) {
                text(x = xa, y = 1, labels = ir_chosen[changes], col = "black")
            }
        }
    }
    dev.off()
}




alpha_col <- function(col, alpha) {
    x <- col2rgb(col) / 255
    return(rgb(x["red", 1], x["green", 1], x["blue", 1], alpha = alpha)    )
}
