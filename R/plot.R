
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
                     lab_size = 2, top_n = 15, wrap_len = 60){
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
    p <- ggplot(datf, aes_string(x = "OR", y = "lp", label = label)) + 
    geom_point(shape = 16, 
               aes_string(size = size, color = color), 
               alpha = alpha) + 
    geom_text_repel(size = lab_size, segment.size = 0.2) + 
    xlab("Odds Ratio") + ylab("-log p") + 
    theme_bw() + 
    guides(color = guide_legend(override.aes = list(size = 1, alpha = 1)))
    return(p)
}

