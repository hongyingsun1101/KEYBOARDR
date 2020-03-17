#' Boundary Table and Decision Matrix for Phase I/II Trials.
#' 
#' Generates boundary table and decision matrix for single-agent phase I/II trials using KEYBOARD design.
#' 
#' @details 
#' The KEYBOARD design relies on the posterior distribution of the toxicity and efficacy to guide dosage transition. To make the decision of dose escalation and
#' de-escalation, given the observed data at the current dose and a boundary table by clinicians, we identify the
#' interval that has the highest joint unit probability mass, which we refer to as
#' \emph{the winner key}. The decision corresponding to the winner key is the decision for dose assignment. 
#' 
#' The prespecified boundary table is determined by the clinicians.  There are four subintervals for toxicity including low, moderate, high, and unacceptable.  The generation of these four intervals require three boundaries provided by the physicians, which can be specified as follows: \cr
#' \enumerate{
#' \item Low interval for toxicity is (0, toxicity.low); toxicity.low the upper boundary for low toxicity interval, which can be specified by using argument \code{toxicity.low}. 
#' \item Moderate interval for toxicity is (toxicity.low, toxicity.moderate); toxicity.low is the upper boundary for low toxicity interval, which can be specified by using argument \code{toxicity.low} and toxicity.moderate is the upper boundary for moderate toxicity interval, which can be specified by using argument \code{toxicity.moderate}.
#' \item High interval for toxicity is (toxicity.moderate, toxicity.high); toxicity.moderate is the upper boundary for moderate toxicity interval, which can be specified by using argument \code{toxicity.moderate} and toxicity.high is the upper boundary for high toxicity interval, which can be specified by using argument \code{toxicity.high}.
#' \item Unacceptable interval for toxicity is (toxicity.high, 1); toxicity.high is the upper boundary for high toxicity interval, which can be specified by using argument \code{toxicity.high}.
#' }
#' 
#' Similarly, there are four subintervals for efficacy as follows:\cr
#' \enumerate{
#' \item Low interval for efficacy is (0, efficacy.low); efficacy.low the upper boundary for low efficacy interval, which can be specified by using argument \code{efficacy.low}. 
#' \item Moderate interval for efficacy is (efficacy.low, efficacy.moderate); efficacy.low is the upper boundary for low efficacy interval, which can be specified by using argument \code{efficacy.low} and efficacy.moderate is the upper boundary for moderate efficacy interval, which can be specified by using argument \code{efficacy.moderate}.
#' \item High interval for efficacy is (efficacy.moderate, efficacy.high); efficacy.moderate is the upper boundary for moderate efficacy interval, which can be specified by using argument \code{efficacy.moderate} and efficacy.high is the upper boundary for high efficacy interval, which can be specified by using argument \code{efficacy.high}.
#' \item Superb interval for efficacy is (efficacy.high, 1); efficacy.high is the upper boundary for high efficacy interval, which can be specified by using argument \code{efficacy.high}.
#' }
#' \cr
#' Suppose that there are d doses in the trial and the current dose is i. Define the number of patients as \eqn{n_i}, the number of patients experiencing toxicity as \eqn{x_i} and the number of responders as \eqn{y_i}.  The trial data can be represented as follows:
#' \deqn{D= {(n_i, x_i, y_i), i = 1, ..., d}} 
#' \cr
#' Assuming the toxicity probability is \eqn{p_i} and the efficacy probability is \eqn{q_i} at dose level i, the probability unit intervals for toxicity and efficacy can be partitioned in to subintervals (a, b) and (c,d). (a, b)  is the subinterval for toxicity probability and (c, d) is the subinterval for efficacy probability. (a, b) x (c,d) is one combination interval.  There are 16 combination intervals in total given 4 subintervals for toxicity rate and 4 subintervals for efficacy intervals.  Clinicians would then be required to provide 16 decisions that correspond to these 16 combination intervals. Decision "D" denotes de-escalation and the next cohort of patients will be treated at the closest lower dose level.  Decision "E" denotes escalation and  the next cohort of patients will be treated at the closest higher dose level.  Decision "S" denotes stay and the next cohort of patients will be treated at the current dose level. The following is an example of a prespecified boundary table with target toxicity rate of 0.2 and target efficacy rate of 0.4:
#' 
#' \tabular{rrrrrrr}{
#'        \tab   \tab Efficacy.low \tab Efficacy.moderate \tab Efficacy.high \tab Efficacy.superb \cr
#'        \tab   \tab  (0,0.25) \tab  (0.25,0.45) \tab (0.45,0.65) \tab (0.65,1) \cr
#'   Toxicity.low \tab  (0,0.15)\tab  E\tab  E\tab  E \tab  E\cr
#'   Toxicity.moderate \tab  (0.15,0.25)\tab  E\tab  E\tab  E \tab  S\cr
#'   Toxicity. high \tab  (0.25,0.35)\tab  D\tab  S\tab  S \tab  S\cr
#'   Toxicity. unaccpetable \tab  (0.35,1.0)\tab  D\tab  D\tab  D \tab D\cr
#' }
#' 
#' For example, the combination interval (0.33, 0.4) x (0,0.2) corresponds to decision "D". This means that the next cohort of patients will be treated at the closest lower level if the observed toxicity rate of current dose falls in (0.33, 0.4)  and the observed efficacy rate falls in (0,0.2). 
#' 
#' Bayesian rule is used to calculate the joint unit probability mass (JUPM) for the toxicity and efficacy combination intervals. For a given combination interval, JUPM is calculated as follows:
#' \deqn{JUPM_{(a,b)}^{(c,d)} =Pr{p_j \in (a,b), q_j \in (c,d) | D} / (b-a)*(d-c) }
#' \eqn{Pr{p_j \in (a,b), q_j \in (c,d) | D} } is the posterior probability of \eqn{p_i} and \eqn{q_i} falling in the subinterval (a,b) and (c,d).  Assume the priors for both \eqn{p_i} and \eqn{q_i} follow independent beta distributions \eqn{beta(\alpha_p,\beta_p)} and \eqn{beta(\alpha_q,\beta_q)} independently. The posterior distributions for \eqn{p_i} and \eqn{q_i} are \eqn{beta(\alpha_p + x_i,\beta_p + n_i - x_i)} and \eqn{beta(\alpha_q + y_i,\beta_q + n_i - y_i)}. Using these posterior distributions, calculate the JUPM for all these 16 combination intervals and find the winning combination interval (a*,b*) and (c*,d*) with the largest JUPM. The decision that corresponds with this winner combination interval is used to treat the next cohort of patients.
#' 
#' 
#' There are two dose exclusion rules applied in the trial design: safety rule and futility rule.
#' Safety rule:
#' if at least 3 patients have been treated at a given dose and given 
#' the observed data the probability of the toxicity rate of
#' the current dose exceeding the target toxicity rate is more
#' than 95\%, we eliminate current and any higher doses from the trial to prevent
#' exposing future patients to unacceptable toxic doses. The probability
#' threshold can be specified with \code{cutoff.eli.toxicity}. When a dose is
#' eliminated, the design recommends the closest lower dose for treating the next cohort of 
#' patients. If the lowest dose is unacceptable toxic, the trial terminates early and
#' no dose is selected as the OBD. This corresponds to a dose assignment of "DUT", which means that de-escalate due to unacceptable high toxicity and exclude the current dose and any dose higher than this dose from the trial.
#' 
#' Futility rule:
#' if at least 3 patients have been treated at a given dose and given
#' the observed data the probability of the efficacy rate of
#' the current dose exceeding the target efficacy rate is less than 30\%, we eliminate this dose from the trial and prevent from
#' exposing future patients to these futile doses. The probability
#' threshold can be specified with \code{cutoff.eli.efficacy}.  This corresponds to two possible dose assignments: "EUE" and "DUE", which both exclude the current dose from the trial. "EUE" denotes escalation due to unacceptable low efficacy and "DUE" denotes de-escalation due to unacceptable low efficacy.
#' 
#' An attractive feature of the KEYBOARD design is that its dose escalation and
#' de-escalation rule can be tabulated before implementing the trial. Thus,
#' when conducting the trial, no real-time calculation or model fitting is needed, and we
#' only need to count the number of patients, the number of DLTs, and the number of responders observed at the current dose and make
#' the decision of dose escalation and de-escalation based on the pre-tabulated
#' decision rules.
#' 
#' @param toxicity.low the upper boundary for low toxicity interval
#' @param toxicity.moderate the upper boundary for moderate toxicity interval
#' @param toxicity.high the upper boundary for high toxicity interval
#' @param efficacy.low the upper boundary for low efficacy interval
#' @param efficacy.moderate the upper boundary for moderate efficacy interval
#' @param efficacy.high  the upper boundary for high efficacy interval
#' @param target.toxicity the target DLT rate
#' @param target.efficacy the target efficacy rate
#' @param cohortsize the cohort size
#' @param ncohort the total number of cohorts
#' @param cutoff.eli.toxicity the cutoff to eliminate an unacceptbale high toxic dose for safety.
#'                            The recommended the default value is (\code{cutoff.eli.toxicity=0.95}).
#' @param cutoff.eli.efficacy the cutoff to eliminate a dose with unacceptable low efficacy.
#'                            The recommended default value is (\code{cutoff.eli.toxicity=0.3})
#' @return \code{get.decision.obd.kb()} returns a prespecified boundary table and a dose assignment decision matrix:
#' \enumerate{
#' \item Prespecified boundary table (\code{$boundary.table})
#' \item Decision matrix (\code{$decision.matrix})
#' 
#' }
#' @author Hongying Sun, Li Tang, Cheng Cheng, and Haitao Pan
#' 
#' @examples
#'  toxicity.low <- 0.15
#'  toxicity.moderate <- 0.25
#'  toxicity.high <- 0.35
#'  efficacy.low <- 0.25
#'  efficacy.moderate <- 0.45
#'  efficacy.high <- 0.65
#'  target.toxicity=0.20
#'  target.efficacy=0.40
#'  cohortsize=3
#'  ncohort=10
#'
#'  decision.obd <- get.decision.obd.kb( toxicity.low = toxicity.low,
#'                  toxicity.moderate= toxicity.moderate,
#'                  toxicity.high = toxicity.high,
#'                  efficacy.low = efficacy.low,
#'                  efficacy.moderate = efficacy.moderate,
#'                  efficacy.high = efficacy.high,
#'                  target.toxicity=target.toxicity,
#'                  target.efficacy=target.efficacy,
#'                  cohortsize=cohortsize, ncohort=ncohort)
#'  print(decision.obd)
#' 
#' @family single-agent phase I/II functions
#' 
#' @note This method is adopted from Li et al (2017)
#' 
#' @references 
#' Daniel H. Li, James B. Whitmore, Wentian Guo and Yuan Ji  Toxicity and Efficacy Probability Interval Design for Phase I Adoptive Cell Therapy Dose-Finding Clinical Trials
#' \emph{Clinical Cancer Research}. 2017; 23:13-20.
#'https://clincancerres.aacrjournals.org/content/23/1/13.long



get.decision.obd.kb <- function(toxicity.low, toxicity.moderate,toxicity.high, efficacy.low, efficacy.moderate,
                                efficacy.high, target.toxicity, target.efficacy, cohortsize, ncohort, cutoff.eli.toxicity= 0.95, cutoff.eli.efficacy=0.3){

    decision.matrix <- matrix(c("E","E","E","E","E","E","E","S","D","S","S","S", "D","D","D","D"),ncol=1, byrow=TRUE)
    efficacy <-matrix(rep( c(0,efficacy.low,efficacy.low,efficacy.moderate,efficacy.moderate,efficacy.high,efficacy.high,1),4), ncol=2, byrow=TRUE)
    toxicity <- matrix( c(0,toxicity.low,0,toxicity.low,0,toxicity.low,0,toxicity.low,toxicity.low,toxicity.moderate,toxicity.low,toxicity.moderate,toxicity.low,toxicity.moderate,toxicity.low,toxicity.moderate,toxicity.moderate,toxicity.high,toxicity.moderate,toxicity.high,toxicity.moderate,toxicity.high,toxicity.moderate,toxicity.high,toxicity.high,1,toxicity.high,1,toxicity.high,1,toxicity.high,1), ncol=2, byrow=TRUE)
    #create decision table
    decision.table <- cbind(toxicity,efficacy, decision.matrix)
    colnames(decision.table) <- c("T1", "T2", "EF1","EF2","DECISION")

    #create output tables to store results
    #to store the total number of patients.
    ntvector <- c()
    ndvector <-  c()
    nrvector <-  c()
    for (i in 1:ncohort) {
        nt<-rep(i*cohortsize,(i*cohortsize+1)^2)
        ntvector <- append(ntvector, nt)
        for (j in 0:(i*cohortsize)){
            nd <- rep(j, (i*cohortsize+1))
            #print(nd)
            ndvector <- append(ndvector, nd)
            nr <- c(0:(i*cohortsize))
            nrvector <- append(nrvector,nr)
        }

    }
    ntmatrix <- matrix(ntvector,ncol=1)
    ndmatrix <- matrix(ndvector,ncol=1)
    nrmatrix <- matrix(nrvector,ncol=1)


    #output table matrix
    outmatrix <- cbind(ntmatrix,ndmatrix,nrmatrix)
    colnames(outmatrix)<- c("N", "T","R")
    outmatrix <- as.data.frame(outmatrix)

    #this function is to return the total number of patients, the number of patients with DLT and the number of patients with reponse given a matrix and an index.
    pinfor <- function(mout, i){
        out <- list(  n = mout[i,1],
                      d = mout[i,2],
                      r = mout[i,3])
        # return (out)
    }
    # given decision table and an index, output the intervals for the toxicity and efficacy and also the dicisions.
    pdecision <- function(decisont, i){
        out <- list (
            t1 = as.numeric(decisont[i,1]),
            t2 = as.numeric(decisont[i,2]),
            e1 = as.numeric(decisont[i,3]),
            e2 = as.numeric(decisont[i,4]),
            d = decisont[i,5]

        )
    }
    # finds the decision that provides the largest jupm
    maxjupm <- function(decision.table, outmatrix, i){
        outputoutm <- pinfor(outmatrix,i)
        nt = outputoutm$n
        dt = outputoutm$d
        rt = outputoutm$r
        jupmv <- c()
        for (rowd in 1:nrow(decision.table)){
            #outdecision <- pdecision(decision.table,rowd)
            a = as.numeric(decision.table[rowd,1])
            b = as.numeric(decision.table[rowd,2])
            c = as.numeric(decision.table[rowd,3])
            d = as.numeric(decision.table[rowd,4])
            #D = outdecision$d
            jupm = round((pbeta(b, (1+dt), (1+nt-dt))-pbeta(a, (1+dt), (1+nt-dt)))*(pbeta(d, (1+rt), (1+nt-rt))-pbeta(c, (1+rt), (1+nt-rt)))/((b-a)*(d-c)),2)
            jupmv <- append(jupmv, jupm)
        }
        jupmv_matrix <- matrix(jupmv, ncol=1)
        colnames(jupmv_matrix)<- "jupm"
        decision.table_jupm <- cbind(decision.table,jupmv_matrix)
        decision.table_jupm <- as.data.frame(decision.table_jupm)
        decision.table_jupm$jupm <- as.numeric(as.character(decision.table_jupm$jupm))
        jupm_max <-max(decision.table_jupm$jupm)
        max_index <- which.max(decision.table_jupm$jupm)
        decision <- decision.table_jupm$D[max_index]

        out <- list(decision.table_jupm=decision.table_jupm, chosen_decision=decision, jupm_max_value =jupm_max)

    }

    #select decision
    decision_picking <-c()
    jupm_picking <-c()
    for (rowoutm in 1:nrow(outmatrix)){
        outjupm <- maxjupm(decision.table, outmatrix, rowoutm)
        decision_p <- outjupm$chosen_decision
        decision_picking <- append(decision_picking,as.character(decision_p))
        jupm_mv <- outjupm$jupm_max_value
        jupm_picking <- append(jupm_picking,jupm_mv)

    }
    decision_matrix <- matrix(decision_picking, ncol=1)
    colnames(decision_matrix)<- "Decision"
    jupm_max_matrix <- matrix(jupm_picking, ncol=1)
    colnames(jupm_max_matrix) <- "Jupm"

    output.matrix <- cbind(outmatrix, decision_matrix,jupm_max_matrix)
    # remove jupm information.
    output.matrix <- output.matrix[,-5]
    output.matrix <- as.data.frame(output.matrix)
    i <- sapply(output.matrix, is.factor)
    output.matrix[i] <- lapply(output.matrix[i], as.character)
    for (rowindex in 1:nrow(output.matrix)){
      if (output.matrix$Decision[rowindex] == "D") {
        if (1- pbeta(target.toxicity,output.matrix$T[rowindex]+1, output.matrix$N[rowindex]-output.matrix$T[rowindex]+1) > cutoff.eli.toxicity) {
          output.matrix$Decision[rowindex] = "DUT"

        } else if (1- pbeta(target.efficacy,output.matrix$R[rowindex]+1, output.matrix$N[rowindex]-output.matrix$R[rowindex]+1)< cutoff.eli.efficacy){
          output.matrix$Decision[rowindex] = "DUE"

        }
      } else if (output.matrix$Decision[rowindex] == "E"){
        if (1- pbeta(target.efficacy,output.matrix$R[rowindex]+1, output.matrix$N[rowindex]-output.matrix$R[rowindex]+1)< cutoff.eli.efficacy){
          output.matrix$Decision[rowindex] = "EUE"
        }

      }
    }



    outlist =list(
            boundary.table=decision.table,
            decision.matrix=output.matrix)

    return (outlist)

}





