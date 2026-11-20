% Run one Dynare estimation of the book's Smets-Wouters model and save
% what export_dynare.m needs. Works in MATLAB and in Octave.
%
% Environment:
%   DYNARE_ROOT  Dynare installation (default C:\dynare\7.0)
%   SW_RUN_DIR   folder holding <SW_FNAME>.mod, the book's three included
%                .mod files and usmodel_data.mat
%   SW_FNAME     Dynare file to run, without extension (sw_table1, sw_long1, ...)
%
% Writes <SW_RUN_DIR>/<SW_FNAME>_results.mat (oo_, M_, options_,
% estim_params_, bayestopt_, elapsed).
root = getenv('DYNARE_ROOT');
if isempty(root), root = 'C:\dynare\7.0'; end
addpath(fullfile(root, 'matlab'));
run_dir = getenv('SW_RUN_DIR');
fname = getenv('SW_FNAME');
cd(run_dir);
if exist('OCTAVE_VERSION', 'builtin'), host = ['Octave ' version]; else, host = ['MATLAB ' version]; end
fprintf('Dynare %s, %s, %s in %s\n', dynare_version(), host, fname, pwd);
t0 = tic;
eval(['dynare ' fname ' noclearall console']);
elapsed = toc(t0);
fprintf('\nWall time: %.1f s\n', elapsed);
out_file = fullfile(run_dir, [fname '_results.mat']);
if exist('OCTAVE_VERSION', 'builtin')
    save('-mat7-binary', out_file, 'oo_', 'M_', 'options_', 'estim_params_', 'bayestopt_', 'elapsed');
else
    save(out_file, 'oo_', 'M_', 'options_', 'estim_params_', 'bayestopt_', 'elapsed');
end
