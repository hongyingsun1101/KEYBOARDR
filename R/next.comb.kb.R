#' Find the Next Dose Combination
#'
#' Determines the dose combination for the next cohort of patients in
#' drug-combination trials that aim to find one MTD.
#'
#' @details
#' Given the observed data thus far, this function determines the dose
#' combination for treating the next cohort of new patients. The observed data
#' are: the number of patients treated at each dose combination (\code{npts}),
#' even if some doses have not been trialed and therefore those numbers are
#' zero; the number of patients who experienced dose-limiting toxicities (DLTs)
#' at each dose combination (\code{ntox}), even if some doses have not been
#' trialed and therefore those numbers are zero; and the level of the current
#' dose (\code{dose.curr}).
#'
#' @param target the target dose-limiting toxicity (DLT) rate.
#' @param npts a \code{J*K} matrix \code{(J<=K)} containing the number of
#'             patients treated at each dose combination.
#' @param ntox a \code{J*K} matrix \code{(J<=K)} containing the number of
#'             patients who experienced a dose-limiting toxicity at each dose
#'             combination.
#' @param dose.curr the current dose combination.
#' @param n.earlystop the early stopping parameter. If the number of patients
#'                    treated at the current dose reaches \code{n.earlystop},
#'                    stop the trial and select the MTD based on the observed
#'                    data.\cr
#'                    The default value, 100, essentially turns off this type
#'                    of early stopping.
#' @param marginL the difference between the target and the left bound of the
#'                "target key" (proper dosing interval) to be defined.\cr
#'                The default is 0.05.
#' @param marginR the difference between the target and the right bound of the
#'                "target key" (proper dosing interval) to be defined.\cr
#'                The default is 0.05.
#' @param cutoff.eli the cutoff to eliminate an overly toxic dose and all
#'                   higher doses for safety.\cr
#'                   The recommended value for general use and default is 0.95.
#' @param extrasafe set \code{extrasafe=TRUE} to impose a more stringent
#'                  stopping rule.\cr
#'                  The default is FALSE.
#' @param offset a small positive number (between 0 and 0.5) to control how
#'               strict the stopping rule is when \code{extrasafe=TRUE}. A
#'               larger value leads to a stricter stopping rule.\cr
#'               The default value of 0.05 generally works well.
#'
#' @return This function returns the recommended dose for treating the next
#'   cohort of patients (\code{$next_dc}).
#'
#' @examples
#' ### Drug-combination trial ###
#'
#' n <- matrix(c(3, 0, 0, 0, 0,
#'               7, 6, 0, 0, 0,
#'               0, 0, 0, 0, 0), ncol=5, byrow=TRUE)
#' y <- matrix(c(0, 0, 0, 0, 0,
#'               1, 1, 0, 0, 0,
#'               0, 0, 0, 0, 0), ncol=5, byrow=TRUE)
#'
#' nxt.comb <- next.comb.kb(target=0.3, npts=n, ntox=y, dose.curr=c(2, 2))
#' summary.kb(nxt.comb)
#'
#' @section Uses:
#' This function uses \code{\link{get.boundary.comb.kb}}.
#'
#' @family drug-combination functions
#'
#' @references
#'
#' 1. Yan F, Mandrekar SJ, Yuan Y. KEYBOARD: A Novel Bayesian Toxicity Probability
#' Interval Design for Phase I Clinical Trials.
#' \emph{Clinical Cancer Research}. 2017; 23:3994-4003.
#' http://clincancerres.aacrjournals.org/content/23/15/3994.full-text.pdf
#' 2. Pan H, Lin R, Yuan Y. Keyboard design for phase I drug-combination trials
#' \emphh{Contemporary Clinical Trials}. 2020
#' https://doi.org/10.1016/j.cct.2020.105972
#'
#'
next.comb.kb <- function(target, npts, ntox, dose.curr, n.earlystop = 100,
                         marginL = 0.05, marginR = 0.05, cutoff.eli = 0.95,
                         extrasafe = FALSE, offset = 0.05) {
    set.seed(1)

    ## simple error checking
    if (npts[dose.curr[1], dose.curr[2]] == 0) {
        stop("Dose entered is not the current dose. \n")
    }
    if (target < 0.05) {
        stop("The target is too low! \n")
    }
    if (target > 0.6) {
        stop("The target is too high! \n")
    }
    if (offset >= 0.5) {
        stop("The offset is too large! \n")
    }
    if (n.earlystop <= 6) {
        stop("Warning: the value of n.earlystop is too low to ensure good operating characteristics. \n",
             "  Recommend n.earlystop = 9 to 18. \n")
    }

    ## obtain dose escalation and de-escalation boundaries
    temp = get.boundary.comb.kb(target, ncohort=150, cohortsize=1, marginL, marginR, cutoff.eli)$boundary
    b.e = temp[2, ]   # escalation boundary
    b.d = temp[3, ]   # de-escalation boundary
    b.elim = temp[4, ]  # elimination boundary

    n = npts
    y = ntox
    earlystop = 0
    d = dose.curr
    nc = n[d[1],d[2]]
    ndose = length(npts)
    elimi = matrix(rep(0, ndose),dim(n)[1],dim(n)[2])  ## indicate whether doses are eliminated

    ## determine if early termination is needed
    if (n[d[1],d[2]] >= n.earlystop) {
        cat("Terminate the trial because the number of patients treated at (", d[1], ", ", d[2], ") has reached",
            n.earlystop,  "\n")
        d = c(99, 99)
        earlystop = 1
    }

    if (!is.na(b.elim[nc])) {
        if (d[1]==1 && d[2]==1 && y[d[1],d[2]]>=b.elim[nc]) {
            d = c(99, 99)
            earlystop = 1
            cat("Terminate the trial because the lowest dose is overly toxic \n")
        }

        ## implement the extra safe rule by decreasing the elimination cutoff for the lowest dose
        if (extrasafe) {
            if (d[1]==1 && d[2]==1 && n[1,1]>=3) {
                if (1-pbeta(target, y[1,1]+1, n[1,1]-y[1,1]+1) > cutoff.eli-offset) {
                    d = c(99, 99)
                    earlystop = 1
                    cat("Terminate the trial because the lowest dose is overly toxic \n")
                }
            }
        }
    }

    ## determine elimination status for combinations
    for (i in 1:dim(n)[1]) {
        for (j in 1:dim(n)[2]) {
            if (n[i,j]>0 && (!is.na(b.elim[n[i,j]]))) {
                if (y[i,j] >= b.elim[n[i,j]]) {
                    elimi[i:dim(n)[1], j:dim(n)[2]] = 1
                }
            }
        }
    }

    out = list("next_dc"=c(NA,NA))
    if (earlystop == 0) {
        ## dose escalation/de-escalation
        if (y[d[1],d[2]] <= b.e[nc]) {
            elevel = matrix(c(1,0,0,1), 2)
            pr_H0 = rep(0, length(elevel)/2)
            nn = pr_H0
            for (i in seq(1, length(elevel)/2, by=1)) {
                ## Code taken from get.oc.comb.kb (for if condition below,
                ## 'n' was originally 'p.true')
                if (d[1] + elevel[1, i] <= dim(n)[1] &&
                    d[2] + elevel[2, i] <= dim(n)[2]) {
                    if (elimi[d[1] + elevel[1, i], d[2] + elevel[2, i]] == 0) {
                        yn = y[d[1] + elevel[1, i], d[2] + elevel[2, i]]
                        nn[i] = n[d[1] + elevel[1, i], d[2] + elevel[2, i]]
                        pr_H0[i] <- pbeta(target+marginR, yn + 0.5, nn[i] - yn + 0.5) -
                                    pbeta(target-marginL, yn + 0.5, nn[i] - yn + 0.5)
                        ## Old code, from BOIN counterpart:
                        #pr_H0[i] <- pbeta(lambda2, yn+0.5, nn[i]-yn+0.5) -
                        #            pbeta(lambda1, yn+0.5, nn[i]-yn+0.5)
                    }
                }
            }
            pr_H0 = pr_H0+nn*0.0005  ## break ties

            if (max(pr_H0)==0) {
                d = d
            }
            else {
                k = which(pr_H0 == max(pr_H0))[as.integer(runif(1)*length(which(pr_H0==max(pr_H0)))+1)]
                d = d+c(elevel[1,k],elevel[2,k])
            }
        }
        else if(y[d[1],d[2]] >= b.d[nc]) {
            delevel = matrix(c(-1,0,0,-1), 2)
            pr_H0 = rep(0, length(delevel)/2)
            nn = pr_H0
            for (i in seq(1, length(delevel)/2, by=1)) {
                if (d[1]+delevel[1,i]>0 && d[2]+delevel[2,i]>0) {
                    yn = y[d[1]+delevel[1,i], d[2]+delevel[2,i]]
                    nn[i] = n[d[1]+delevel[1,i], d[2]+delevel[2,i]]
                    ## Code taken from get.oc.comb.kb:
                    pr_H0[i] = pbeta(target + marginR, yn + 0.5, nn[i] - yn + 0.5) -
                               pbeta(target - marginL, yn + 0.5, nn[i] - yn + 0.5)
                    ## Old code, from BOIN counterpart:
                    # pr_H0[i] = pbeta(lambda2, yn+0.5, nn[i]-yn+0.5) -
                    #            pbeta(lambda1, yn+0.5, nn[i]-yn+0.5)
                }
            }
            pr_H0 = pr_H0+nn*0.0005 ## break ties

            if (max(pr_H0)==0) {
                d = d
            }
            else {
                k = which(pr_H0==max(pr_H0))[as.integer(runif(1)*length(which(pr_H0==max(pr_H0)))+1)]
                d = d+c(delevel[1,k],delevel[2,k])
            }
        }
        else {
            d = d
        }
        out = list("next_dc"=d)
    }
    return(out)
}
