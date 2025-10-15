#' Collect Nodes
#'
#' @description
#' Collect results run on parallel nodes
#'
#' @details If `N` simulations were run across `k` nodes, then each node will have
#' the results from approximately `N/k` simulations. This function gathers the the locations
#' of all results into a single list. When you start to run your simulations in parallel, this function
#' works best if the directory does not already contain results from a previous attempt to run your
#' simulations. Therefore, if you ran your `N` simulations and it has failed, remove the `k`
#' subdirectories and start again.
#'
#' @param clusters The number of nodes used in parallel.
#' @param directory The directory where the `k` subdirectories are found. Each `k` subdirectory
#' is named after its node.
#'
#' @return List of `k` filenames with full paths.
#' @export
#'
#' @examples
#' # load libraries for parallelization
#' library(parallel)
#' library(pbapply)
#' # load inputs from data
#' data(inputs)
#' # Example of input list to generate a two-arm-two-stage trial with binary outcome data
#' # Run 4 simulations with lapply without parallelization
#' num.sims <- 4
#' conf <- lapply(
#' 1:num.sims,
#' runSingleTrial,
#' inputs=inputs,
#' save.plot=FALSE,
#' print=FALSE,
#' directory = '')
#' # Now run in parallel across two nodes
#' clusters <- 2
#' cl <- makeCluster(clusters)
#' directory <- tempdir()
#' res.list <- pblapply(
#' 1:num.sims,
#' runSingleTrial,
#' inputs=inputs,
#' save.plot=FALSE,
#' directory=directory,
#' cl=cl)
#' stopCluster(cl)
#' # now collect the results from the two nodes
#' res.list <- collectNodes(clusters, directory)
#' res.list <- do.call("rbind", res.list)
#' res.list <- data.frame(res.list)
#' # now your results are in one list

collectNodes <- function(clusters, directory){
  # exclude parent
  dirs = list.dirs(directory)[-1]
  # get basenames
  basenames = lapply(dirs, basename)
  # find basenames that are numeric only
  numfiles = suppressWarnings(lapply(basenames, as.numeric))
  use.dirs = 1 - is.na(numfiles)
  # mask list of directories
  dirs = dirs[which(!is.na(numfiles))]
  # get the time stamp
  times=lapply(dirs,function(x){
    info = file.info(x)
    t.str = strptime(info$ctime, "%Y-%m-%d %H:%M:%S")
    round(as.numeric(format(t.str, "%H")) +
            as.numeric(format(t.str, "%M"))/60, 2)
  } )
  # find which ones were created at the same time, in the cluster
  dirs = dirs[which(times == unique(times)[unlist(lapply(unique(times), function(x) {
    sum(times==x)==clusters
  }))][[1]])]
  # gather all csvs from dirs
  csvs = unlist(lapply(dirs, function(x){dir(x, full.names=T, pattern=".csv+") }))
  return(lapply(csvs, read.csv))
}
