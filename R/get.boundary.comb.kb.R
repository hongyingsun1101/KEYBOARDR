#' Dose Escalation and De-escalation Boundaries for Drug-combination Trials
#'
#' Generates the optimal dose escalation and de-escalation boundaries for
#' conducting a drug-combination trial with the KEYBOARD design.
#'
#' @details
#' The KEYBOARD design relies on the posterior distribution of the toxicity
#' probability to guide dosage. To make the decision of dose escalation and
#' de-escalation, given the observed data at the current dose, we identify the
#' interval that has the highest posterior probability, which we refer to as
#' the "strongest key". This key represents where the true dose-limiting
#' toxicity (DLT) rate of the current dose is most likely located. If the
#' strongest key is located on the left side of the "target key", we escalate
#' the dose (because it means that the observed data suggests that the current
#' dose is most likely to represent under-dosing); if the strongest key is
#' located on the right side of the target key, we de-escalate the dose
#' (because the data suggests that the current dose represents overdosing); and
#' if the strongest key is the target key, we retain the current dose (because
#' the observed data supports that the current dose is most likely to be in the
#' proper dosing interval).
#' Graphically, the strongest key is the one with the largest area under the
#' posterior distribution curve of the DLT rate of the current dose.
#'
#' \figure{keyboard.jpg}
#' 
#' An attractive feature of the KEYBOARD design is that its dose escalation and
#' de-escalation rule can be tabulated before the onset of the trial. Thus,
#' when conducting the trial, no calculation or model fitting is needed, and we
#' only need to count the number of DLTs observed at the current dose and make
#' the decision of dose escalation and de-escalation based on the pre-tabulated
#' decision rules.
#'
#' Given all observed data, we use matrix isotonic regression to obtain the
#' estimate of the toxicity rate of the combination of dose level j of drug A
#' and dose level k of drug B, and select the MTD as the combination with the
#' toxicity estimate that is closest to the target. When there are ties, we
#' randomly choose one as the MTD.
#'
#' For patient safety, we apply the following Bayesian overdose control rule
#' after each cohort:
#' if at least 3 patients have been treated at the given dose and
#' the observed data indicate that the probability of the toxicity rate of
#' the current combination dose being above the target toxicity rate is more
#' than 95%, we eliminate that and higher doses from the trial to prevent
#' exposing future patients to these overly toxic doses. The probability
#' threshold can be specified with \code{cutoff.eli}. If the lowest dose
#' combination (1, 1) is overly toxic, the trial terminates early and no dose
#' is selected as the MTD.
#'
#' @param target The target dose-limiting toxicity (DLT) rate.
#' @param ncohort A scalar specifying the total number of cohorts in the trial.
#' @param cohortsize The number of patients in the cohort.
#' @param marginL The difference between the target and the left bound of the
#'                "target key" (proper dosing interval) to be defined.\cr
#'                The default is 0.05.
#' @param marginR The difference between the target and the right bound of the
#'                "target key" (proper dosing interval) to be defined.\cr
#'                The default is 0.05.
#' @param cutoff.eli The cutoff to eliminate an overly toxic dose and all
#'                   higher doses for safety.\cr
#'                   The recommended value for general use and default is 0.95.
#'
#' @return The function returns a matrix, which includes the dose escalation
#'   and de-escalation boundaries, as well as the elimination boundary.
#'
#' @note In most clinical applications, the target DLT rate is often a rough
#'   guess, but finding a dose level with a DLT rate reasonably close to the
#'   target rate (which ideally would be the MTD) is what interests the
#'   investigator.
#'
#' @examples
#' ### Drug-combination trial ###
#'
#' bound <- get.boundary.comb.kb(target=0.3, ncohort=10, cohortsize=3)
#' print(bound)
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
get.boundary.comb.kb <- function(target, ncohort, cohortsize,n.earlystop=100,
                                 marginL=0.05, marginR=0.05, cutoff.eli=0.95,offset=0.05,extrasafe=TRUE) {
  ## Get cutoffs for keys
  getkeys <- function(a1, a2)
  {
    delta=a2-a1
    lkey=NULL; rkey=NULL
    i=0
    cutoff=0.3
    while (cutoff>0)
    {
      i=i+1
      cutoff = a1-i*delta
      lkey = c(cutoff, lkey)
    }
    lkey[lkey<0]=0

    i=0; cutoff=0.3
    while (cutoff<1)
    {
      i=i+1
      cutoff = a2+i*delta
      rkey = c(rkey, cutoff)
    }
    rkey[rkey>1]=1
    key=c(lkey, a1, a2, rkey)

    return(key)
  }

  ## Identify the key with the highest posterior tocutoff.elicity probability
  keys = getkeys(target - marginL, target + marginR)
  nkeys = length(keys) - 1

  npts = ncohort * cohortsize
  targetp = rep(NULL, nkeys)

  a = b = 1  # Hyperparameters used in beta prior
  decision_table = matrix(NA, nrow = npts + 1, ncol = npts)
  for (ntr in 1:npts) {
    elim = 0
    for (ntox in 0:ntr) {

        if (ntox >= 3 & 1 - pbeta(target, ntox + a, ntr - ntox + b) > cutoff.eli) {
          elim = 1
          break
        }

      # targetp = rep(NULL, nkeys)
      for(i in 1:nkeys) {
        comp = 1
        if (i == 1 || i == nkeys) {
          # Compensation factor for incompleted keys:
          comp = (marginL + marginR) / (keys[i+1] - keys[i])
        }
        targetp[i] = (pbeta(keys[i+1], ntox + a, ntr - ntox + b) -
                      pbeta(keys[i], ntox + a, ntr - ntox + b)) * comp
      }
      highkey = max(which(targetp == max(targetp)))
      targetkey = which(keys == (target - marginL))

      if (highkey > targetkey) {
        decision_table[ntox+1, ntr] <- "D"
      }
      if (highkey == targetkey) {
        decision_table[ntox+1, ntr] <- "S"
      }
      if (highkey < targetkey) {
        decision_table[ntox+1, ntr] <- "E"
      }
    }
    if (elim == 1) {
      decision_table[(ntox+1):(ntr+1), ntr] <- rep("DU", ntr - ntox + 1)
    }
  }
  colnames(decision_table) <- 1:npts
  rownames(decision_table) <- 0:npts

  boundary = matrix(NA, nrow = 4, ncol = npts)
  boundary[1, ] = 1:npts
  for (i in 1:npts) {
    if (length(which(decision_table[, i] == "E"))) {
      boundary[2, i] = max(which(decision_table[, i] == "E")) - 1
    }
    else {
      boundary[2, i] = -1
    }
    if (length(which(decision_table[, i] == "D"))) {
      boundary[3, i] = min(which(decision_table[, i] == "D")) - 1
    }
    else if (length(which(decision_table[, i] == "DU"))) {
      boundary[3, i] = min(which(decision_table[, i] == "DU")) - 1
    }
    if (length(which(decision_table[, i] == "DU"))) {
      boundary[4, i] = min(which(decision_table[, i] == "DU")) - 1
    }
  }
  colnames(boundary) <- c(rep("", npts))
  rownames(boundary) <- c("Number of patients treated",
                          "Escalate if # of DLT <=",
                          "de-escalate if # of DLT >=",
                          "Eliminate if # of DLT >=")


  if (extrasafe) {
    stopbd = NULL
    ntrt = NULL
    for (n in 1:npts) {
      ntrt = c(ntrt, n)
      if (n < 3) {
        stopbd = c(stopbd, "NA")
      }
      else {
        for (ntox in 1:n) {
          if (1 - pbeta(target, ntox + 1, n - ntox +
                        1) > cutoff.eli - offset) {
            stopneed = 1
            break
          }
        }
        if (stopneed == 1) {
          stopbd = c(stopbd, ntox)
        }
        else {
          stopbd = c(stopbd, "NA")
        }
      }
    }
    stopboundary = data.frame(rbind(ntrt, stopbd)[, 1:min(npts, n.earlystop)])
    rownames(stopboundary) = c("The number of patients treated at the lowest dose  ",
                               "Stop the trial if # of DLT >=        ")
    stopboundary = data.frame(stopboundary)
    #colnames(stopboundary) = rep("", min(npts, n.earlystop))
    colnames(stopboundary) = 1:dim(stopboundary)[2]

  }



 if(extrasafe){
   return(list(boundary=boundary,safe=stopboundary))
 }else{
   return(list(boundary=boundary))
 }


}

###  obtain dose escalation/de-escalation boundaries up to n=20 for conducting the trial
get.boundary.comb.kb(target=0.2, ncohort=8, cohortsize=3, marginL=0.05, marginR=0.05,
                    cutoff.eli=0.95, offset=0.05,extrasafe=TRUE)
get.boundary.comb.kb(target=0.2, ncohort=8, cohortsize=3, marginL=0.05, marginR=0.05,
                     cutoff.eli=0.95, offset=0.05,extrasafe=FALSE)



#   target=0.2; ncohort=8; cohortsize=3; marginL=0.05; marginR=0.05; cutoff.eli=0.95; offset=0.05;extrasafe=TRUE
