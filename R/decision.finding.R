#' Returns dose escalation and de-escalation given any time point of the trial
#'
#' Outputs the decision given the decision table, the number of total patients, the number of patients experiencing toxicity and the number of responders
#'
#' @param out.matrix a decision matrix
#' @param n the number of enrolled subjects
#' @param t the number of subjects experienced dose limiting toxicities (DLT)
#' @param r the number of responders
#'
#' @return \code{decision.finding()} returns the decision.
#'
#' @examples
#'  n=6
#'  t=1
#'  r=3
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
#'  output.matrix <- get.decision.obd.kb( toxicity.low = toxicity.low,
#'                  toxicity.moderate= toxicity.moderate,
#'                  toxicity.high = toxicity.high,
#'                  efficacy.low = efficacy.low,
#'                  efficacy.moderate = efficacy.moderate,
#'                  efficacy.high = efficacy.high,
#'                  target.toxicity=target.toxicity,
#'                  target.efficacy=target.efficacy,
#'                  cohortsize=cohortsize, ncohort=ncohort)$decision.matrix
#' decision <- decision.finding (out.matrix=output.matrix, n=n, t=t, r=r)
#' decision
#'

decision.finding <- function(out.matrix, n, t, r){
  rowindex <- which(out.matrix$N==n & out.matrix$T==t & out.matrix$R == r)
  decision <- out.matrix$Decision[rowindex]
  decision <- as.character(decision)
  return (decision)
}
