// Long MH chain for the book's Smets-Wouters model (one of several run in
// parallel). Copy to sw_long<k>.mod, set the seed, and put it in a folder
// with the book's sw_defs.mod, sw_model_def.mod, sw_estimation_def.mod,
// usmodel_data.mat and the mode file of the book's run
// (sw_table1/Output/sw_table1_mode.mat).
@#include "sw_defs.mod"
@#include "sw_model_def.mod"
@#include "sw_estimation_def.mod"

set_dynare_seed(2007);   // a different seed for each chain (2007 + k used)

// The book's estimation options, started at the book run's posterior mode
// and Hessian; one block per process, so the chain starts at the mode.
estimation(
	nograph,
	nodiagnostic,
	datafile=usmodel_data,
	mode_compute=0,
	mode_file=sw_table1_mode,
	first_obs=71,
	presample=4,
	lik_init=2,
	prefilter=0,
	mh_replic=150000,
	mh_nblocks=1,
	mh_jscale=0.20,
	mh_drop=0.3333333333
);
