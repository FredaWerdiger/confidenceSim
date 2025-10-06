# Functions used in the design and simulation of confidence trials


#' Get Group Sequential Design
#' @description Generate boundaries for a group sequential design using a one-sided test
#' @details To generate confidence-based thresholds, we are interested in a one-sided test and alpha is set at 0.025.
#' To generate the stopping thresholds, specify either looks or information rates, and the alpha spending function if
#' difference from O-Brien-Fleming-type.
#'
#' @param info.rates Analysis times expressed as rate of information accrual. Expects a vector
#' with the last item representing the final analysis and equal to 1.
#' For example, information rates for two-stage trial with interim analysis half way through is c(0.5, 1).
#' One of two options for expressing analysis times. Either 'info.rates' or 'looks' must be specified.
#' @param looks Analysis times expressed by number of patients accrued at each point. Expects a vector
#' with the last item being equal to the maximum sample size. For example.
#' looks for a three stage trial with maximum sample size of 300 and analysis planned every 100 patients is c(100, 100, 100).
#' One of two options for expressing analysis times. Either 'info.rates' or 'looks' must be specified.
#' @param as.type Time of alpha spending function to use. Options are as outlined by 'rpact'.
#' Default is 'asOF' \(O'Brien-Fleming-type\).
#'
#' @return Returns an 'rpact' TrialDesign object.
#' @seealso [rpact::getDesignGroupSequential()]
#' @export
#'
#' @examples
#' # calculate critical values for a two-stage trial with an interim analysis half-way through
#' # Use Pocock-type alpha spending
#' design <- getGSDesign(info.rates = c(0.5, 1), as.type = 'asP')
#' critical.stagewise.alpha.levels <- design$stageLevels
#'
#' # calculate values for a 6-stage trial with a maximum sample size of 1000
#' # interim analysis begins at 500 patients accrued and continues every 100 patients after
#' design <- getGSDesign(looks = seq(500, 1000, 100))
#'

getGSDesign <- function(info.rates=NULL, looks=NULL,  as.type="asOF"){

  if (is.null(looks) & is.null(info.rates)){
    stop("Looks or information rates is required.")
  }
  if (!is.null(looks)) {
    num.looks = length(looks)
    t = looks/max(looks) #
  }
  if (!is.null(info.rates)){
    t = info.rates
  }
  if (! as.type %in% c("OF", "P", "WT", "PT", "HP", "WToptimum", "asP", "asOF", "asKD",
                       "asHSD", "asUser", "noEarlyEfficacy")){
    print(pastt0("Alpha spending option ",as.type, " not available. Changing to default (asOF)."))
    as.type = "asOF"
  }
  # using a one-sided test for confidence threshold
  alpha = 0.025
  sided = 1
  design = rpact::getDesignGroupSequential(
    sided=sided,
    alpha=alpha,
    informationRates = t,
    typeOfDesign = as.type
  )
  return(design)
}


#' Get Confidence Levels from Group Sequential Bounds
#' @description Derive confidence-based decision thresholds for efficacy and inferiority from a group sequential design.
#' @param design TrialDesign object generated from 'getGSDesign'.
#' @seealso [getGSDesign()]
#'
#' @return List of values:
#' \itemize{
#' \item{critical.values}{Critical z-levels at each stage}
#' \item{alpha.spending.one.sided}{One-sided alpha critical values for each stage}
#' \item{alpha.spending.cumulative}{One-sided alpha accumulated by each stage}
#' \item{confidence.threshold.efficacy}{Upper bounds for confidence i.e. declare efficacy if exceeded}
#' \item{confidence.threshold.inferiority}{Lower bounds for confidence i.e. declare inferiority if below}
#' }
#'
#' @export
#'
#' @examples
#' # confidence bounds for a 6-stage trial with a maximum sample size of 1000
#' bounds <- getConfidenceFromBounds(getGSDesign(looks=seq(500, 1000, 100)))
#'

getConfidenceFromBounds <- function(design){
  if (typeof(design)!="environment"){
    stop("Function expecting TrialDesign object.")
  }
  # for efficacy
  confidence.threshold = 1 - design$stageLevels
  # for inferiority, confidence thresholds are p-values
  bounds=list(critical.values=design$criticalValues,
              alpha.spending.one.sided=design$stageLevels,
              alpha.spending.cumulative=design$alphaSpent,
              confidence.threshold.efficacy=confidence.threshold,
              confidence.threshold.inferiority=design$stageLevels)
  return(bounds)
}


#' Get Critical Bounds from Confidence Thresholds
#'
#' @description
#' Derive traditional frequentist critical values from frequentist confidence thresholds for confidence in treatment benefit.
#'
#' @details
#' During a confidence trial, efficacy and inferiority is determined by the level of confidence in treatment benefit.
#' Efficacy is declared if this confidence level exceeds a pre-specified boundary,
#' and inferiority is declared if this confidence levels falls below a second pre-specified valye.
#' Given confidence-based thresholds for efficacy and inferiority, and the
#' sidedness of the test, this function returns the traditional frequentist p-value.
#'
#'
#' @param num.treat.arms Number of treatment arms (excludes control). Default is 2.
#' @param conf.lower Confidence in treatment benefit boundary for inferiority i.e. stop for inferiority if confidence in benefit is below this.
#' Default is 0.01.
#' @param conf.upper Confidence in treatment benefit boundary for efficacy i.e. stop for efficacy if confidence in benefit is above this.
#' Default is 0.99.
#' @param p.sided Sidedness of statistical test, 1 (one-sided) and 2 (two-sided). Default is 1.
#'
#' @return List of values:
#' \itemize{
#' \item{conf.lower: confidence in treatment benefit lower bound}
#' \item{z.score.lower: critical value corresponding to lower confidence bound}
#' \item{p.value.lower: p-value correponding to lower confidence bound}
#' \item{conf.upper: confidence in treatment benefit upper bound}
#' \item{z.score.upper: critical value corresponding to upper confidence bound}
#' \item{p.value.upper: p-value corresponding to upper confidence bound}
#' \item{p.value: sidedness of test}
#' }
#'
#' @export
#'
#' @examples
#' # Running the function on default values
#' bounds <- getBoundsFromConfidence()
#'
#' # to make adjustments for multiple arms
#' bounds <- getBoundsFromConfidence(num.treat.arms = 3)
#'
getBoundsFromConfidence <- function(num.treat.arms=2,
                                    conf.lower=0.01,
                                    conf.upper=0.99,
                                    p.sided=1
){
  # If there are more than three treatment arms, the lower threshold is lower for each arm
  if(num.treat.arms <=2){
    conf.lower = conf.lower
  } else if (num.treat.arms >2){
    conf.lower = conf.lower/2
  }
  critical.value.upper = qnorm(conf.upper)
  critical.value.lower = qnorm(conf.lower)
  if(p.sided == 1){
    p.value = 'one-tailed'
    p.upper = 1 - conf.upper
    p.lower = conf.lower
  } else if(p.sided == 2){
    p.value = 'two-tailed'
    p.upper = (1 - conf.upper)*2
    p.lower = conf.lower * 2
  }
  bounds = list(conf.lower=conf.lower, z.score.lower=critical.value.lower,
                p.value.lower=p.lower, conf.upper=conf.upper, z.score.upper=critical.value.upper,
                p.value.upper=p.upper, p.value=p.value)
  return(bounds)
}


# get a parameter list to generate trial
#' Get Parameter List
#' @description
#' Generate a parameter list to generate a frequentist confidence trial
#'
#' @param looks Vector of analysis times expressed by either number of patients accrued at each point or by rate of information accumulated.
#' If the former, last item should be the maximum sample size.
#' If the latter, last item is 1.
#' Expects a vector with length equal to the total number of total looks.
#' @param nmax Maximum sample size, specified if information rates are used for 'looks'.
#' @param perpetual Whether to run the trial perpetually (TRUE) or not (FALSE).
#' If TRUE, new treatment arms will be added when treatment arms are dropped, until there are no more arms left.
#' All treatments are included via the 'resprate' parameter. Default is FALSE.
#' @param alloc.ratio Allocation ratios for study arms relative to each other.
#' Expects vector with length equal to number of arms including control.
#' First number corresponds to control ratio.
#' @param num.per.block Number from each arm per block, for blocked randomization to balance co-variates.
#' Block size is 'sum(num.per.block)'. If a single number is provided, it will be assume to apply to each arm.
#' @param final.visit The number of days after intervention when the response information becomes available.
#' @param as.type The type of alpha spending function to use in group sequential design.
#'  Default is 'asOF', O'Brien-Fleming-type.
#' @param alpha The alpha threshold to apply to each pairwise comparison to control in the final analysis.
#' Used together with the 'MONITOR FUTILITY', when an alpha spending function is not needed.
#' Default is 0.05, assuming a two-sided test.
#' @param multiarm.mode For multiple treatment arms, describes how arms are evaluated at each stage:
#' \itemize{
#' \item{"CONFIDENCE-BASED"(default): Evaluate arms against confidence-based rules }
#' \item{"DROP WORST": Drop the worst performing arm, and carry the remaining promising arms}
#' \item{"SELECT BEST": Select the best performing arm to carry forward, drop the rest}
#' \item{"ALL PROMISING": Carry forward all promising arms}
#' \item{"MONITOR FUTILITY": Only monitor for futility}
#' }
#' @param lmb
#' @param lmb.conf.thresh
#' @param outcome.type
#' @param estimator.type
#' @param resprate
#' @param ppm
#' @param special
#'
#' @return
#' @export
#'
#' @examples
getparlist = function(looks=seq(500,1000,100),
                      nmax=NULL,
                      perpetual=FALSE,
                      alloc.ratio=c(1,1),
                      num.per.block=c(1,1),
                      final.visit=0,
                      as.type="asOF",
                      alpha=0.05,
                      multiarm.mode="CONFIDENCE-BASED",
                      lmb=0.10,
                      lmb.conf.thresh=0.9,
                      outcome.type='BINARY',
                      estimator.type='odds ratio',
                      resprate=c(0.3,0.5),
                      ppm=rep(15, 300),
                      special=NULL) {
  # design parameters
  # looks: a vector interim analysis according to data accrued; max(looks) is the maximum sample size for a treatment
  # perpetual: whether or not the trial continutes beyond nmax to test more arms (can arms be added?)
  # ppm: vector with number of patients per month. Length of vector is number of months.
  # alloc.ratio: vector for allocation ratios, control is the first number
  # num.per.block: for blocked randomisation to balance covariates
  # outcome.type: Options are BINARY, CONTINUOUS, ORDINAL and TRIPLETS
  # resprate: The response rate for control and then each treatment arm.
  #           # for binary outcome a vector of probabilities is expected
  # for continuous outcome, a vector of (mean, sd) is expected FOR EACH arm
  # for ordinal outcome only taking interact rates
  # can take win-loss-tie triplets for generating log-odd
  # There can be more resprates than arms allocated. This mean arms are added once others are dropped.
  # final visit: in days, the time for the response information to become available
  #lmb: Define lack of meaningful benefit, relative to whatever your outcome metric is
  # as.type: the type of alpha-spending function to use. for RPACT.
  # multiarm.mode: for multiple treatment arms, how to evaluated at each stage
  #                 options: "DROP WORST" "SELECT BEST" "ALL PROMISING" "MONITOR FUTILITY"
  # Special = passs something special


  # check if looks are given as time rate
  if (max(looks) == 1){
    if (!is.null(nmax)){
      looks = looks * nmax
    } else{
      stop("Maximum sample size not supplied.")
    }
  }

  nmax = max(looks)

  parlist = list(nmax = nmax, looks=looks, num.looks = length(looks))

  if (!is.null(special)){
    parlist$misc=special
  } else {parlist$misc='None'}

  parlist$alloc.ratio = alloc.ratio
  if (length(num.per.block) != length(alloc.ratio)){
    num.per.block = rep(num.per.block[1], length(alloc.ratio))
  }
  parlist$num.per.block = num.per.block

  parlist$final.visit = final.visit # default assumes immediate follow-up

  if (!toupper(multiarm.mode) %in% c("ALL PROMISING", "SELECT BEST", "DROP WORST", "CONFIDENCE-BASED", "MONITOR FUTILITY")){
    print(paste0("Mulitarm mode option ", multiarm.mode, " not available. Changing to default (CONFIDENCE-BASED)"))
    multiarm.mode = "CONFIDENCE-BASED"
  }

  parlist$multiarm.mode = toupper(multiarm.mode)
  parlist$alpha=alpha

  # if only monitoring for futility
  if (parlist$multiarm.mode == 'MONITOR FUTILITY'){
    parlist$as.type = NULL
    parlist$critical.values = NULL
    parlist$alpha.spending = NULL
    parlist$confidence.bounds.efficacy = NULL
    parlist$confidence.bounds.inferiority = NULL
  } else {

    # decision thresholds via rpact
    if (! as.type %in% c("OF", "P", "WT", "PT", "HP", "WToptimum", "asP", "asOF", "asKD",
                         "asHSD", "asUser", "noEarlyEfficacy")){
      print(pastt0("Alpha spending option ",as.type, " not available. Changing to default (asOF)."))
      as.type = "asOF"
    }

    parlist$as.type = as.type
    design = getGSDesign(looks=looks, as.type=as.type)
    bounds = getConfidenceFromBounds(design)
    parlist$critical.values = bounds$critical.values
    parlist$alpha.spending = bounds$alpha.spending
    parlist$confidence.bounds.efficacy = bounds$confidence.threshold.efficacy
    parlist$confidence.bounds.inferiority = bounds$confidence.threshold.inferiority
  }

  # LMB // Lack of Meaningful Benefit
  parlist$lmb.threshold = lmb
  parlist$lmb.confidence.threshold = lmb.conf.thresh

  # response rates and accrual
  if (!toupper(outcome.type) %in% c("BINARY", "ORDINAL", "CONTINUOUS", "TRIPLETS")){
    print(paste0("Outcome type ", outcome.type, " not in options."))
  }
  parlist$outcome.type = toupper(outcome.type)

  # check if there are enough response rates
  if (length(resprate) < length(alloc.ratio)){
    stop("Not enough response rates to allocated arm.")
  }

  if (is.numeric(resprate)){
    # if the input is a list of numbers, then it is binary
    parlist$resprate = resprate
    parlist$outcome.type = toupper(outcome.type)
    if(outcome.type != "BINARY"){
      print("Outcome type changed to BINARY according to numeric response rate")
      parlist$outcome.type = "BINARY"
    }
  } else if (is.character(resprate)){
    # built in interact response rates
    if(outcome.type != "ORDINAL"){
      print("Outcome type changed to ORDINAL according to string response rate")
      parlist$outcome.type = "ORDINAL"
    }
    if (resprate == 'interact'){
      parlist$resprate = strokeTrials::getInteractResprate()
    } else if (resprate == 'interact.rev'){
      parlist$resprate = strokeTrials::getOppositeResprate()
    } else if (resprate == 'interact.same'){
      parlist$resprate = strokeTrials::getSameResprate()
    }
  } else if (is.list(resprate)){
    # check if each has two items (mean and std)
    if (length(resprate[1])==2){
      if (outcome.type != "CONTINUOUS") {
        print("Outcome type changed to CONTINUOUS according to list of mean/std response rate")
        parlist$outcome.type = "CONTINUOUS"}
    } else {
      if(outcome.type != "ORDINAL"){
        if (outcome.type != "TRIPLETS"){
          print("Outcome type changed to ORDINAL according to list response rate")
          parlist$outcome.type = "ORDINAL"}
      }
    }
    parlist$resprate = resprate
  }

  # check binary has estimator type
  if (parlist$outcome.type == 'BINARY'){
    if (grepl('risk', tolower(estimator.type)) | grepl('diff', tolower(estimator.type))){
      parlist$estimator.type = 'risk diff'
    } else {
      parlist$estimator.type = 'odds ratio'
    }
  }

  # accrual
  parlist$ppm = ppm
  parlist$perpetual=perpetual

  # check if there are more response rates than allocations
  if ((length(resprate) == length(alloc.ratio)) & perpetual==TRUE){
    print("No extra treaments provided. Turning off perpertual.")
    parlist$perpetual = FALSE
  }

  resp.str =  sapply(seq(length(parlist$resprate)), function(x){
    paste(unlist(sapply(parlist$resprate[[x]]*100, round,0)), collapse='-')})
  parlist$resprate.str = resp.str
  parlist$outfilename = paste0(
    length(parlist$alloc.ratio), "arm",
    length(looks), "stage",
    parlist$outcome.type,
    "lmb-", 100* parlist$lmb.threshold,
    parlist$as.type, "-", parlist$misc)
  return(parlist)
}


# generate patient accrual with Poisson distribution
# how many patients do we accrue to get the number of follow-ups?
# is this trial running perpetually?
getAccrual = function (numsubjects, ppm, follow.up=0,
                       cont.recruit=FALSE,
                       perpetual=FALSE){

  if (!(is.numeric(numsubjects) & is.numeric(ppm))){
    ppm = suppressWarnings(as.numeric(ppm))
    numsubjects = suppressWarnings(as.numeric(numsubjects))
    if (is.na(ppm) | is.na(numsubjects)){
      stop("Inputs numsubjects or ppm not numeric.")
    }
  }

  monthin = NULL
  ptsin = 0  # we care about this if cont.recruit = FALSE
  ptsout = 0 # we care about this is cont.recruit = TRUE
  for (i in 1:length(ppm)){
    # cumsum means that each item in the vector accumulates on the previous entry
    # this helps them represent TIME continuously
    temp = cumsum(rexp(n=((numsubjects*2)-ptsin), rate=ppm[i])) # sample extra patients
    temp = temp[temp <= 1] # only this months patients
    # combine patients in this months with total patients
    monthin = c(monthin, (i - 1) + temp)
    ptsin <- length(monthin)
    ptsout <- length(monthin[monthin > follow.up]) # will equal pts in if follow.up = 0
    if (!perpetual){
      if(ptsout >= numsubjects){
        # do we halt recruiting
        if (cont.recruit==FALSE){
          monthin <- monthin[1:numsubjects]
        } else {
          # do we continue recruiting after reaching max
          monthin <- monthin[monthin <= monthin[numsubjects] + follow.up]
        }
        break
      }
    }
  }
  return(monthin) # get a vector size of number of patients you need to get numsubjects followup
}

# function to randomize patients
# Treatment 1: 1 to one block, one to another
# each block gets an equal distribution of treatments to remove the effect
# that is different between blocks.
getBlockedArm = function(numsubjects, num.per.block, prob=NULL){

  if (!is.numeric(numsubjects)){
    numsubjects = suppressWarnings(as.numeric(numsubjects))
    if (is.na(numsubjects)){
      stop("'numsubjects' not numeric.")
    }
  }

  num.per.block = as.vector(num.per.block)
  num.arms = length(num.per.block)

  if (num.arms <= 1 | is.na(num.arms) ){stop('Check "num.per.block".')}

  block.size = sum(num.per.block)
  # if no probabilities are supplied, make all arms equal
  if (is.null(prob)){
    prob=rep(1, num.arms)
  }
  # if some arms probabilities are set to zero change block size and need
  # TODO: size
  if (any(prob==0)){
    block.size =sum(num.per.block[prob>0])
  }
  need = ceiling(numsubjects/block.size)
  # making sure there is a particular distribution
  as.vector(replicate(
    need,sample(rep.int(1:num.arms, num.per.block),
                prob=rep.int(prob, num.per.block),
                size=block.size)))[1:numsubjects]
}
####################e
# GENERATE OUTCOMES
###################

# binary outcome expects a probability response rate for each arm
# it should be supplied thus c(control rate, treatment1 rate,..., treatmentn rate)
# 1 is generate at that specified rate.
getDataBin = function(arm, resprate){
  if (is.na(resprate[arm])){
    stop("No response rate for this treatment.")
  }
  rbinom(n=1, size=1, prob=resprate[arm])}

# ordinal outcomes (e.g. mRS) expect a list of probabilities for each arm
getDataOrd <- function(arm, resprate){
  if (is.na(resprate[arm])){
    stop("No response rate for this treatment.")
  }
  which(( rmultinom(
    n=1,
    size=1,
    prob = resprate[[arm]]))==1)
}

# continuous outcomes expect a mean and standard deviation for each arm
# to generate values from a distribution (e.g., normal)
getDataCont <- function(arm, resprate, dist='norm'){
  if (is.na(resprate[arm])){
    stop("No response rate for this treatment.")
  }
  if (dist=='norm'){
    rnorm(1,
          mean=resprate[[arm]][1],
          sd=resprate[[arm]][2])
  }

}



# function to data that is available at interim
# Set as.followup = TRUE is you want n to be the amount of available follow-up
getCurrentData = function(datlist, looktime, n, as.followup=TRUE){
  n = suppressWarnings(as.numeric(n))
  if (is.na(n)){
    stop("Num subjects (n) supplied as  NA.")
  }
  if (is.na(looktime)){
    stop("Looktime supplied as NA.")
  }

  if (as.followup){
    # the time of interest is when the 500th reaches follow up
    looktime = datlist$obstime[n]
    if (is.na(looktime)){
      stop(paste0("There is no ", n, "th subject in the dataset"))
    }
    enrolled = (datlist$arrival.day <= looktime)
  } else{
    # the time of interest is when the 500th is enrolled
    if (max(datlist$arrival.day) < looktime){
      stop("Latest arrival day is less than chosen looktime.")
    }

    enrolled = (datlist$arrival.day <= looktime) & (datlist$subjid <= n)
  }
  # create subset of data based on bool
  # x[y] is the function where x is the original list and y is the boolean array
  newdatlist = lapply(datlist, FUN=function(x,y) x[y], y=enrolled)
  # data where we have the follow-up/outcome data
  # if there is no time to follow up then all of them are known
  newdatlist$known = newdatlist$obstime<=looktime
  return(c(newdatlist, looktime=looktime))
}

## Get sufficient statistics
getSuffStats = function(datlist) {
  num.enrolled = table(datlist$arm)
  num.known = with(datlist, tapply(datlist$known, datlist$arm, sum))
  num.unknown = num.enrolled - num.known
  remarms = sort(unique(datlist$arm))
  if(length(unique(datlist$dat)) < 3) {
    num.resp = tapply(datlist$dat[datlist$known], datlist$arm[datlist$known], sum)
    num.fail = num.known - num.resp
    resprate = num.resp/num.known
    suffstats = list(
      num.enrolled = num.enrolled,
      num.known = num.known,
      num.unknown = num.unknown,
      num.resp = num.resp,
      num.fail = num.fail,
      resprate = resprate)
    suffstats$formattedrate = paste0(num.resp, "/", num.known,
                                     "(", resprate, ")")
  } else if (length(unique(datlist$dat)) < 30){
    # ordinal data
    num.resp.all = sapply(remarms, function(x) {
      resp = datlist$dat[datlist$known][datlist$arm[datlist$known]==x]
      do.call("rbind",
              lapply(1:length(unique(datlist$dat)),
                     function(y){length(which(resp==y))}))
    }, USE.NAMES = FALSE)
    cols = c("control",paste0("treatment", seq(ncol(num.resp.all)-1)))
    resprate = asplit(prop.table(num.resp.all ,2), MARGIN=2)
    num.resp = asplit(num.resp.all, MARGIN=2)
    names(num.resp) = cols
    names(resprate) = cols
    suffstats = list(
      num.enrolled = num.enrolled,
      num.known = num.known,
      num.unknown = num.unknown,
      num.resp = num.resp,
      resprate = resprate)
    formattedrate = lapply(1:length(resprate), function(x){sapply(resprate[[x]]*100, round,0)})
    names(formattedrate) = cols
    # TODO: format the following for print
    suffstats$formattedrate = formattedrate

  }else {
    means = sapply(remarms, function(x) {
      mean(datlist$dat[datlist$known][datlist$arm[datlist$known]==x])})
    sds =  sapply(remarms, function(x) {
      sd(datlist$dat[datlist$known][datlist$arm[datlist$known]==x])})

    suffstats=list(
      mean = means,
      sd = sds,
      arms = remarms,
      num.enrolled = num.enrolled,
      num.known = num.known,
      num.unknown = num.unknown
    )
    suffstats$formattedrate=paste0(signif(means, 2))
  }

  return(suffstats)
}
