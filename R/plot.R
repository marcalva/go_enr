
#' Plot enrichments with bubbles
#'
#' Takes data frame output from fet_go or fet_react, and plots 
#' the log fold change against the p-value. The point size 
#' corresonds to the number of terms
#' 
#' @param datf Data frame returned by fet_go or fet_react.
#' @param p_thresh Only terms with adjusted p-value less than this 
#'  value are plotted.
#' @param size Column label for point size.
#' @param label Label for the point. Can be none or any column in data frame.
#' @param color Column label to color points.
#' @param wrap_len Wrap the labels so every line is at most this many 
#'  characters long.
#'
#' @return ggplot object
plot_enr <- function(datf, p_thresh = 0.05, size = "GenesInTerm", 
                     label = "Name", color = NULL, alpha = 0.2, 
                     lab_size = 2, top_n = 15, wrap_len = 60, 
                     add_vline = FALSE, vline_col = "red"){
    require(ggplot2)
    require(ggrepel)
    k <- datf[,"p_adj"] < p_thresh
    if (sum(k) == 0){
        warning("No singificant terms")
    }
    datf <- datf[k,,drop = FALSE]
    datf[,label] <- as.character(datf[,label])
    if (sum(k) > 0){
        datf[,label] <- sapply(datf[,label], function(s) { paste(strwrap(s, wrap_len), collapse="\n") } )
    }
    if (top_n < nrow(datf)){
        nr <- min(top_n, nrow(datf))
        datf[(nr+1):nrow(datf), label] <- ""
    }
    datf[,"lp"] <- -log10(datf[,"p_adj"])
    p <- ggplot(datf, aes_string(x = "OR", y = "lp", label = label))
    if (add_vline) p <- p + geom_vline(xintercept = 0, color = "red", linetype="dashed")
    p <- p + geom_point(shape = 16, 
               aes_string(size = size, color = color), 
               alpha = alpha) + 
    geom_text_repel(size = lab_size, segment.size = 0.2) + 
    xlab("Odds Ratio") + ylab("-log p") + 
    theme_bw() + 
    guides(color = guide_legend(override.aes = list(size = 1, alpha = 1)))
    return(p)
}

# reverse and log transform a vector of values
log10_rev_trans <- function(x){
    trans_new("log10_rev",
              function(x) {x = log10(rev(x)); x[is.na(x)] = log10(0.05); return(x)},
              function(x) {x = rev(10 ^ (x)); x[is.na(x)] = 0.05; return(x)},
              log_breaks(10),
              domain = c(1e-100, Inf))
}

log10_rev_breaks <- function(x, n=5){

    rng <- range(x, na.rm = TRUE)
    lxmin <- floor(log10(rng[1])) + log10(2)
    lxmax <- ceiling(log10(rng[2])) - log10(2)

    lby <- floor((lxmax - lxmin)/n) + 1

    breaks <- rev(10^seq(lxmax, lxmin, by=-lby))
    # breaks <- 10^seq(lxmin, lxmax, lby)
    return(breaks)
} 
    
format_power10 <- function(x){
    x <- signif(x, digits = 2)
    lapply(x, function(y){
           pow_num <- floor(log10(y))
           base_num <- y * 10^-pow_num
           # if (pow_num <= 2 && pow_num >= -2) bquote(.(y))
           # else bquote(.(base_num) %*% 10^.(pow_num))
           bquote(.(base_num) %*% 10^.(pow_num))
              })
}

log10_rev_trans <- function(x){
    trans <- function(x) -log(x, 10)
    inv <- function(x) 10^(-x)
    trans_new("log10_rev", trans, inv, breaks = log10_rev_breaks, domain = c(1e-100, Inf))
    trans_new("log10_rev", trans, inv, breaks = log10_rev_breaks,
              format = format_power10, domain = c(1e-100, Inf))
    # format = function(x) math_format(10^.x),
}

heatmap_enr <- function(lor_df, p_df, max_p = 0.05){
    require(ggplot2)
    require(scales)
    require(reshape2)

    # melt data frames
    lor_df <- lor_df[rownames(p_df), colnames(p_df)]
    p_dfm <- reshape2::melt(as.matrix(p_df))
    lor_dfm <- reshape2::melt(as.matrix(lor_df))

    lor_dfm[,"p"] <- p_dfm[,"value"]
    # colnames(lor_dfm)[colnames(lor_dfm) == "value"] <- "LOR"

    # remove entries with p > max_p
    kp <- lor_dfm[,"p"] <= max_p
    lor_dfm <- lor_dfm[kp,,drop=FALSE]
    
    minp <- min(lor_dfm[,"p"], na.rm = TRUE)
    maxv <- max(abs(lor_dfm[,"value"]), na.rm = TRUE)

    p <- ggplot(lor_dfm, aes_string(x = "Var1", y = "Var2", color = "value", size = "p")) + 
        geom_point() + 
        theme_bw() + 
        scale_size(trans="log10_rev", limits = c(0.05, minp)) + 
        scale_color_distiller(palette = "RdYlBu",
                              name = "log\nodds\nratio",
                              limits = c(-1, 1) * maxv) + 
        theme(text = element_text(size = 14),
              axis.title = element_blank(),
              axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
              axis.ticks =  element_blank(),
              panel.border = element_blank()) +
        guides(size = guide_legend(override.aes = list(shape = 1, color = "black")))

    return(p)
}

