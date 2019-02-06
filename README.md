This folder contains all necessary files to fit the binomial regression model with spatio-temporal uncertainty described in Heaton et al (2019).  The contents of this folder are as follows:
	-- FitBronchiolitisData.R: function for fitting the model
	-- AMCMCUpdate.R: function that updates Metropolis proposal variances
	-- CensusData.R: census information for each grid box in Norfolk VA
	-- GetAdjMatrix.R: function for defining adjacency matrices
	-- Norfolk200mgrid.Rdata: grid centroids for Norfolk VA at 200m level
	-- NorfolkFullPredictionsPM25.Rdata: predictions and SE's of PM 2.5 at grid points in Norfolk
	-- NorfolkFullPredictionsSO2.Rdata: predictions and SE's of SO2 at grid points in Norfolk
	-- SimulatedCaseInfo.txt: a simulated dataset of bronchiolitis cases in Norfolk VA
	-- SimulatedControlInfo.txt: a simulated dataset of bronchiolitis controls in Norfolk VA
	-- README.md: README file containing instructions on how to run the code

FitBronchiolitisData.R is the main function for fitting the model.  This function uses the following libraries:
	-- LatticeKrig version 7.0
	-- data.table version 1.11.2
	-- parallel version 3.5.1
To run the code, the user will need to specify the following variables (all specified at the top of the file):
	-- t.jit.len - the temporal jittering length (in number of days)
	-- grid.box.len - the 1/2 size of the discretization grid of the spatial domain
	-- s.jit.len - the spatial jittering distance
	-- include.X.unc - if TRUE the model will include uncertainty in the pollution measurements in model fitting
	-- include.ST.unc - if TRUE the model will include spatio-temporal uncertainty in model fitting
	-- ncores.unc.update - if include.ST.unc=TRUE, the number of cores used to perform the sampling of true spatio-temporal locations
	-- ncores.param.update - number of cores used to calculate the likelihood of model parameters
Parallel processing is done via mclapply() from the parallel library.  Note that the mclapply() function is incompatible with Intel MKL libraries hence users will need to set MKL_NUM_THREADS=1 in the command line to take advantage of parallel processing.

FitBronchiolitisData.R requires approximately 20 minutes to prep the data by calculating overlap probabilities with the spatial grid cells (actual time depends on the size of the jittering lengths and number of observations) and calculate basis functions via eigendecompositions. After the initial spin up period, each iteration of the MCMC algorithm takes, approximately, 7 seconds using only a single core but gets as low as 2 seconds per iteration using 20 cores.