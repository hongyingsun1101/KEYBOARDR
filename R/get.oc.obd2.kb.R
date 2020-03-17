#'
#' Generate Operating Characteristics for OBD Finding
#'
#' Report operating characteristics for optimal biological dose (OBD) finding
#'
#' @param target.toxicity the target toxicity rate
#' @param target.efficacy the target efficacy rate
#' @param ncohort the total number of cohorts
#' @param cohortsize the cohortsize
#' @param n.early the early stopping parameter. If the number of patients treated at
#'                the current dose reaches \code{n.early}, stop the trial
#'                and select the optimal biological dose (OBD) based on the observed data. The default
#'                value \code{n.early=100} essentially turns off the type
#'                of early stopping.
#' @param startdose the starting dose level
#' @param p.true a vector containing the true toxicity probabilities of the
#'              investigational dose levels.
#' @param q.true a vector containing the true efficacy  probabilities of the
#'              investigational dose levels.
#' @param ntrial the total number of trials to be simulated
#' @param seed the random seed for simulation
#' @param p1 the cutoff lower limit for safety utility function 1, described in the Details section.
#' @param p2 the cutoff upper limit for safety utility function 1, described in the Details section.
#' @param q1 the cutoff lower limit for efficacy utility function 1, described in the Details section.
#' @param q2 the cutoff upper limit for efficacy utility function 1, described in the Details section.
#' @param cutoff.eli.toxicity the cutoff to eliminate a dose with unacceptable high toxicity for safety.
#'                           The recommended value is (\code{cutoff.eli.toxicity=0.95}).
#' @param cutoff.eli.efficacy the cutoff for the futility rule, the acceptable lowest efficacy.
#'                            The recommended value is  (\code{cutoff.eli.toxicity=0.3})
#' @param w1.toxicity the weight for toxicity utility function 2 and 3, described in the Details section.
#' @param w2.toxicity the weight for toxicity utility function 3, described in the Details section.
#' @param indicator the indicator cutoff for utility function 3, described in the Details section.
#'
#' @details A large number of  trials are simulated to characterize the operating characteristics of the KEYBOARD design under the prespecified true toxicity probabilities and true efficacy probabilities of the investigational doses. Dose assignment rule follows the rule described in the function \code{get.decision.obd.kb()}.
#'
#' There are stopping rules the KEYBOARD design:
#' \enumerate{
#' \item stop the trial if the lowest dose is eliminated from the trial due to high unacceptable toxicity.
#' \item stop the trial if the number of patients treated at current dose is larger than or equal to \code{n.earlystop}.
#' }
#'
#' @return \code{get.oc.obd2.kb()} returns the operating characteristics of the KEYBOARD design as a list, including:\cr
#' \enumerate{
#'    \item selection percentage at each dose level using utility function 1 (\code{$selpercent1}), \cr
#'    \item selection percentage at each dose level using utility function 2 (\code{$selpercent2}), \cr
#'    \item selection percentage at each dose level using utility function 3 (\code{$selpercent3}), \cr
#'    \item the number of patients treated at each dose level (\code{$npatients}), \cr
#'    \item the number of dose limiting toxicities (DLTs) observed at each dose level (\code{$ntox}), \cr
#'    \item the number of responders observed at each dose level (\code{$neff}), \cr
#'    \item the average number of DLTs (\code{$totaltox}), \cr
#'    \item the average number of responders (\code{$totaleff}), \cr
#'    \item the average number of patients (\code{$totaln}), \cr
#'    \item the percentage of early stopping without selecting the OBD using utility function 1 (\code{$percentstop1}), \cr
#'    \item the percentage of early stopping without selecting the OBD using utility function 2 (\code{$percentstop2}), \cr
#'    \item the percentage of early stopping without selecting the OBD using utility function 3 (\code{$percentstop3}), \cr
#'    \item data.frame (\code{$simu.setup}) containing simulation parameters, such as target, p.true, etc.
#'
#'
#' }
#'
#'
#' @examples
#' target.toxicity <- 0.30
#' target.efficacy <- 0.40
#' p.true <-c(0.08,0.30,0.60,0.80)
#' q.true <- c(0.25,0.40,0.25,0.50)
#' oc.obd2.kb <- get.oc.obd2.kb(target.toxicity=target.toxicity,
#'               target.efficacy= target.efficacy, ncohort=20,
#'               cohortsize= 3,  p.true= p.true, q.true= q.true)
#' oc.obd2.kb
#' summary.kb(oc.obd2.kb)
#' plot.kb(oc.obd2.kb)
#' plot.kb(oc.obd2.kb$selpercent1)
#' plot.kb(oc.obd2.kb$selpercent2)
#' plot.kb(oc.obd2.kb$selpercent3)
#' plot.kb(oc.obd2.kb$npatients)
#' plot.kb(oc.obd2.kb$ntox)
#' plot.kb(oc.obd2.kb$neff)
#'
#' @family single-agent phase I/II functions
#'
#' @references
#' Daniel H. Li, James B. Whitmore, Wentian Guo and Yuan Ji  Toxicity and Efficacy Probability Interval Design for Phase I Adoptive Cell Therapy Dose-Finding Clinical Trials
#' \emph{Clinical Cancer Research}. 2017; 23:13-20.
#'https://clincancerres.aacrjournals.org/content/23/1/13.long
#'
#' Suyu Liu, Valen E. Johnson. A robust Bayesian dose-finding design for phase I/II clinical trials
#' \emph{Biostatistics}. 2016; 17:249-263.
#' https://academic.oup.com/biostatistics/article/17/2/249/1744018
#'
#' Yanhong Zhou, J. Jack Lee, Ying Yuan. A utility-based Bayesian optimal interval (U-BOIN) phase I/II design to identify the optimal biological dose for targeted and immune therapies
#' \emph{Statistics in Medicine}. 2019; 38:S5299-5316.
#' https://onlinelibrary.wiley.com/doi/epdf/10.1002/sim.8361
#'
get.oc.obd2.kb <- function( target.toxicity, target.efficacy,ncohort=10, cohortsize=3, n.early=100,
                           startdose=1, p.true, q.true, ntrial = 1000, seed = 6, p1=0.15, p2=0.40, q1=0.3, q2=0.6,cutoff.eli.toxicity= 0.95,
                           cutoff.eli.efficacy=0.3, w1.toxicity =0.33, w2.toxicity=1.09, indicator =  target.toxicity){

    # Error checking
    if(target.toxicity<0.05) {cat("Error: the target toxcicity rate is too low! \n"); return();}
    if(target.toxicity>0.6)  {cat("Error: the target toxicity rate is too high! \n"); return();}
    if(target.efficacy<= 0.19) {cat("Error: the target efficacy rate is too low! \n");  return();}



    set.seed(seed)
    ndose = length(p.true);
    npts = ncohort * cohortsize;
    #toxicity outcome matrix
    Y=matrix(rep(0, ndose*ntrial), ncol=ndose);
    #efficacy outcome matrix
    E=matrix(rep(0, ndose*ntrial), ncol=ndose);
    # matrix to store the total number of patients
    N=matrix(rep(0, ndose*ntrial), ncol=ndose);
    # matrix to store selected dose level
    dselect1 = rep(0, ntrial)
    dselect2 = rep(0, ntrial)
    dselect3 = rep(0, ntrial)
    dearlystop = rep(0, ntrial)

    # get decision table
    decision.matrix.output <- get.decision.obd2.kb(target.toxicity=target.toxicity, target.efficacy=target.efficacy, cohortsize=cohortsize, ncohort=ncohort)$decision.matrix

    # this function outputs the decision given the dose-finding table, the number of total patients, the number of patients experiencing toxicity and the number of responders
    decision.finding <- function(out.matrix, n, t, r){
        rowindex <- which(out.matrix$N==n & out.matrix$T==t & out.matrix$R == r)
        decision <- out.matrix$Decision[rowindex]
        decision <- as.character(decision)
        return (decision)
    }


    # simulation trials
    for (trial in 1:ntrial){
        # the number of patients experiencing toxity at each level
        y <- rep(0, ndose);
        # the number of patients experiencing reponse at each level
        e <- rep(0, ndose);
        # the number of total patients treated at each level
        n <- rep(0, ndose);
        # an indicator to check whehter the trial terminates early
        earlystop = 0;
        d=startdose;
        # an indicator to check whehter the doses are eliminated
        elimi = rep(0, ndose);

        for (i in 1:ncohort){
            y[d] = y[d] + sum(runif(cohortsize) < p.true[d]);
            e[d] = e[d] + sum(runif(cohortsize) < q.true[d]);
            n[d] = n[d] + cohortsize;
            if(n[d] >= n.early) break;
            #get dose decision matrix
            decision.result = decision.finding(decision.matrix.output, n[d], y[d], e[d])
            # early stop rules
            if (n[d] >=3 ){
                if (1- pbeta(target.toxicity, y[d]+1, n[d]-y[d]+1) > cutoff.eli.toxicity){
                    elimi[d:ndose]=1;
                    if( d==1) {earlystop=1; break;}
                }
            }
            if (n[d] >=3 ){
                if (1- pbeta(target.efficacy, e[d]+1, n[d]-e[d]+1) < cutoff.eli.efficacy){
                    elimi[d]=1;
                }
            }
            if(!is.null(decision.result)){

                # dose transition
                if(decision.result == "E" && d != ndose && elimi[d+1]==0 )  {
                    d=d+1;
                } else if(decision.result == "D" && d != 1 && elimi[d-1]==0) {
                    d=d-1;
                } else if (decision.result == "S") {
                    d = d ;
                } else if (decision.result == "DUT" && d != 1 && elimi[d-1]==0){
                    d=d-1;
                    elimi[d:ndose]=1;
                } else if (decision.result=="DUE" && d != 1 && elimi[d-1]==0){
                    d=d-1;
                    elimi[d]=1;
                } else if (decision.result == "EUE" && d != ndose && elimi[d+1]==0){
                    d=d+1;
                    elimi[d]=1;
                } else { d=d; }
            }



        }
        Y[trial,] = y;
        E[trial,] = e;
        N[trial,] = n;
        if(earlystop ==1){
            dselect1[trial]=99;
            dselect2[trial]=99;
            dselect3[trial]=99;
            dearlystop[trial]=9999;

        } else{
            dselect1[trial]= select.obd.kb( target.toxicity =target.toxicity, target.efficacy=target.efficacy, npts = n, ntox = y, neff = e)$obd1
            dselect2[trial]= select.obd.kb( target.toxicity =target.toxicity, target.efficacy=target.efficacy, npts = n, ntox = y, neff = e)$obd2
            dselect3[trial]= select.obd.kb( target.toxicity =target.toxicity, target.efficacy=target.efficacy, npts = n, ntox = y, neff = e)$obd3
        } # use select.obd function.
    }

    # output results
    selpercent1 = rep(0, ndose);
    selpercent2 = rep(0, ndose);
    selpercent3 = rep(0, ndose);
    nptsdose = apply(N,2, mean);
    ntoxdose = apply(Y, 2, mean);
    neffdose = apply(E,2, mean);
    for (i in 1:ndose){
        selpercent1[i] = sum(dselect1==i)/ntrial*100;
        selpercent2[i] = sum(dselect2==i)/ntrial*100;
        selpercent3[i] = sum(dselect3==i)/ntrial*100;
    }
    out=list(name = "get.oc.obd.kb",  ## to identify object for summary.kb() function.
             selpercent1=selpercent1,
             selpercent2=selpercent2,
             selpercent3=selpercent3,
             npatients=nptsdose,
             ntox=ntoxdose,
             neff=neffdose,
             totaltox=sum(Y)/ntrial,
             totaleff=sum(E)/ntrial,
             totaln=sum(N)/ntrial,
             earlystop = sum(dearlystop==9999)/ntrial*100,
             percentstop1=sum(dselect1== 99)/ntrial*100,
             percentstop2=sum(dselect2== 99)/ntrial*100,
             percentstop3=sum(dselect3== 99)/ntrial*100,
             simu.setup=data.frame(target.toxicity=target.toxicity, target.efficacy=target.efficacy,p.true=p.true, q.true =q.true,
                                   ncohort=ncohort, cohortsize = cohortsize,
                                   startdose = startdose,
                                   ntrial = ntrial, dose=1:ndose)
    );
    return (out)
}
