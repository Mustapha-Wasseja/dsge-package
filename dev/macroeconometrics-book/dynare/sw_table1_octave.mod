// The book's sw_table1.mod for Octave: identical except mode_compute.
// The book uses mode_compute=1 (MATLAB's fmincon with
// optim=('Algorithm','active-set')), which Octave does not have;
// mode_compute=4 (csminwel) is Dynare's default. Save as sw_table1.mod in
// a folder with the book's three included .mod files and usmodel_data.mat.
@#include "sw_defs.mod"
@#include "sw_model_def.mod"
@#include "sw_estimation_def.mod"

estimation(
	nograph,
	nodiagnostic,
	datafile=usmodel_data,
	mode_compute=4,
	first_obs=71,
	presample=4,
	lik_init=2,
	prefilter=0,
	mh_replic=2000,
	mh_nblocks=2,
	mh_jscale=0.20,
	mh_drop=0.2
);
