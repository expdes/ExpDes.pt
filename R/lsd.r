#' Comparacao multipla: LSD
#'
#' \code{lsd} Realiza o teste t (LSD) para comparacao multipla
#' de medias.
#' @param y Vetor numerico ou complexo contendo a variavel
#' resposta.
#' @param trt Vetor numerico ou complexo contendo os
#' tratamentos.
#' @param DFerror Grau de liberdade do residuo.
#' @param SSerror Soma de quadrados do residuo.
#' @param alpha Significancia do teste de Bootstrap.
#' @param group TRUE ou FALSE
#' @param main Titulo
#' @return E retornada a comparacao das medias segundo o teste
#' LSD.
#' @author Eric B Ferreira,
#'  \email{eric.ferreira@@unifal-mg.edu.br}
#' @author Denismar Alves Nogueira
#' @author Portya Piscitelli Cavalcanti
#' @seealso \code{\link{snk}}, \code{\link{duncan}},
#' \code{\link{ccboot}}, \code{\link{lsdb}},
#' \code{\link{scottknott}}, \code{\link{tukey}},
#' \code{\link{ccF}}.
#' @examples
#' data(ex1)
#' attach(ex1)
#' dic(trat, ig, quali = TRUE, mcomp = "lsd", sigT = 0.05)
#' @importFrom "stats" "AIC" "coef" "dnorm" "fitted"
#' "fitted.values" "nls" "ppoints" "pt" "qnorm" "qt"
#' "residuals"
#' @export

lsd<-function (y, trt, DFerror, SSerror, alpha = 0.05, group = TRUE,    main = NULL)
{
    MSerror <- SSerror/DFerror
    name.y <- paste(deparse(substitute(y)))
    name.t <- paste(deparse(substitute(trt)))
    junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
    means <- tapply.stat(junto[, 1], junto[, 2], stat = "mean")
    sds <- tapply.stat(junto[, 1], junto[, 2], stat = "sd")
    nn <- tapply.stat(junto[, 1], junto[, 2], stat = "length")
    means <- data.frame(means, std.err = sds[, 2]/sqrt(nn[, 2]),
        replication = nn[, 2])
    names(means)[1:2] <- c(name.t, name.y)
    ntr <- nrow(means)
    #Tprob <- qtukey(1 - alpha, ntr, DFerror)
    Tprob <- qt(1-(alpha/2),DFerror)*sqrt(2)
    nr <- unique(nn[, 2])
    nfila <- c("Alpha", "Error Degrees of Freedom", "Error Mean Square",
        "Critical Value of Studentized Range")
    nvalor <- c(alpha, DFerror, MSerror, Tprob)
    #cat("\nStudy:", main)
    #cat("\n\nHSD Test for", name.y, "\n")
    xtabla <- data.frame(...... = nvalor)
    row.names(xtabla) <- nfila
    #print(xtabla)
    #cat("\nTreatment Means\n")
    #print(data.frame(row.names = NULL, means))
    if (group) {
        if (length(nr) == 1) {
            HSD <- Tprob * sqrt(MSerror/nr)
            #cat("\nHonestly Significant Difference", HSD)
        }
        else {
            nr1 <- 1/mean(1/nn[, 2])
            HSD <- Tprob * sqrt(MSerror/nr1)
            #cat("\nHonestly Significant Difference", HSD)
            #cat("\nHarmonic Mean of Cell Sizes ", nr1)
            #cat("\n\nDifferent HSD for each comparison")
        }
        cat("\nTeste t (LSD)\n------------------------------------------------------------------------")
        cat("\nGrupos  Tratamentos  Medias\n")
        output <- order.group(means[, 1], means[, 2], means[,4], MSerror, Tprob, means[, 3], parameter = 0.5)
        cat('------------------------------------------------------------------------\n')
    }
    if (!group) {
        comb <- combn(ntr, 2)
        nn <- ncol(comb)
        dif <- rep(0, nn)
        pvalue <- rep(0, nn)
        for (k in 1:nn) {
            i <- comb[1, k]
            j <- comb[2, k]
            dif[k] <- abs(means[i, 2] - means[j, 2])
            sdtdif <- sqrt(MSerror * (1/means[i, 4] + 1/means[j,
                4]))
            pvalue[k] <- round(1 - ptukey(dif[k] * sqrt(2)/sdtdif,
                ntr, DFerror), 4)
        }
        tr.i <- comb[1, ]
        tr.j <- comb[2, ]
        #cat("\nComparison between treatments means\n\n")
        #print(data.frame(row.names = NULL, tr.i, tr.j, diff = dif, pvalue = pvalue))
        output <- data.frame(trt = means[, 1], means = means[,2], M = "", N = means[, 4], std.err = means[, 3])
    }
#    return(output)
}
