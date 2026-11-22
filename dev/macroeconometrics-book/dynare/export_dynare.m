% Export a finished Dynare estimation of the book's Smets-Wouters model
% to CSV, and evaluate Dynare's log-likelihood and log-prior at every
% parameter vector listed in points.csv (if present).
%
% Environment:
%   SW_RUN_DIR   folder holding <fname>.mod and its Dynare output
%   SW_FNAME     Dynare file name without extension (sw_table1, sw_long, ...)
%   SW_RESULTS   .mat file with oo_, M_, options_, estim_params_, bayestopt_
%   SW_OUT       output folder for the CSV files
%   SW_THIN      keep every SW_THIN-th MH draw in <fname>_draws.csv (default 1)
%   DYNARE_ROOT  Dynare installation (default C:\dynare\7.0)
% Works in MATLAB and in Octave.
root = getenv('DYNARE_ROOT');
if isempty(root), root = 'C:\dynare\7.0'; end
addpath(fullfile(root, 'matlab'));
dynare_config();                 % puts estimation/ etc. on the path
run_dir = getenv('SW_RUN_DIR');
fname   = getenv('SW_FNAME');
out     = getenv('SW_OUT');
cd(run_dir);
load(getenv('SW_RESULTS'));      % oo_ M_ options_ estim_params_ bayestopt_
if ~exist(out, 'dir'), mkdir(out); end

names = bayestopt_.name(:);      % Dynare's order: stderr first, then params
nx = numel(names);
is_sd = false(nx, 1);
is_sd(1:estim_params_.nvx) = true;
names(is_sd) = regexprep(names(is_sd), '^stderr\s+', '');   % 'stderr eta_a' -> 'eta_a'
key = names;
key(is_sd) = strcat('sd_e.', names(is_sd));

% ---- posterior mode and the Hessian-based standard deviations ----------
mode_x = zeros(nx, 1); mode_sd = nan(nx, 1);
for i = 1:nx
    if is_sd(i)
        mode_x(i) = oo_.posterior_mode.shocks_std.(names{i});
        mode_sd(i) = oo_.posterior_std_at_mode.shocks_std.(names{i});
    else
        mode_x(i) = oo_.posterior_mode.parameters.(names{i});
        mode_sd(i) = oo_.posterior_std_at_mode.parameters.(names{i});
    end
end

% ---- MH summaries (if an MH run exists) --------------------------------
has_mh = isfield(oo_, 'posterior_mean');
stat = nan(nx, 5);               % mean, median, sd, hpd_inf, hpd_sup
if has_mh
    flds = {'posterior_mean', 'posterior_median', 'posterior_std', ...
            'posterior_hpdinf', 'posterior_hpdsup'};
    for i = 1:nx
        grp = 'parameters'; if is_sd(i), grp = 'shocks_std'; end
        for j = 1:numel(flds)
            f = flds{j};
            if isfield(oo_, f) && isfield(oo_.(f), grp) && isfield(oo_.(f).(grp), names{i})
                stat(i, j) = oo_.(f).(grp).(names{i});
            end
        end
    end
end

% bounds from estimated_params (only stderr and plain parameters here)
assert(estim_params_.nvn == 0 && estim_params_.ncx == 0 && estim_params_.ncn == 0);
lb = [estim_params_.var_exo(:, 3); estim_params_.param_vals(:, 3)];
ub = [estim_params_.var_exo(:, 4); estim_params_.param_vals(:, 4)];
fid = fopen(fullfile(out, [fname '_summary.csv']), 'w');
fprintf(fid, 'parameter,dynare_name,lower,upper,mode,sd_at_mode,mean,median,sd,hpd_inf,hpd_sup\n');
for i = 1:nx
    fprintf(fid, '%s,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n', ...
            key{i}, names{i}, lb(i), ub(i), mode_x(i), mode_sd(i), stat(i, :));
end
fclose(fid);

% ---- scalars ------------------------------------------------------------
fid = fopen(fullfile(out, [fname '_scalars.csv']), 'w');
fprintf(fid, 'quantity,value\n');
if isfield(oo_, 'MarginalDensity')
    md = oo_.MarginalDensity;
    if isfield(md, 'LaplaceApproximation'), fprintf(fid, 'log_mdd_laplace,%.17g\n', md.LaplaceApproximation); end
    if isfield(md, 'ModifiedHarmonicMean'), fprintf(fid, 'log_mdd_mhm,%.17g\n', md.ModifiedHarmonicMean); end
end
if isfield(oo_, 'posterior') && isfield(oo_.posterior, 'optimization') && ...
        isfield(oo_.posterior.optimization, 'log_density') && ~isempty(oo_.posterior.optimization.log_density)
    fprintf(fid, 'logpost_at_mode,%.17g\n', oo_.posterior.optimization.log_density);
end
if exist('elapsed', 'var'), fprintf(fid, 'wall_seconds,%.17g\n', elapsed); end
opt_names = {'mh_replic', 'mh_nblck', 'mh_jscale', 'mh_drop', 'first_obs', ...
             'presample', 'lik_init', 'nobs'};
for j = 1:numel(opt_names)
    if isfield(options_, opt_names{j}) && isnumeric(options_.(opt_names{j})) && ...
            isscalar(options_.(opt_names{j}))
        fprintf(fid, '%s,%.17g\n', opt_names{j}, options_.(opt_names{j}));
    end
end
fclose(fid);

% ---- raw MH draws, every block ------------------------------------------
mhdir = fullfile(run_dir, fname, 'metropolis');
if has_mh && exist(mhdir, 'dir')
    rec = load(fullfile(mhdir, [fname '_mh_history_0.mat']));
    fid = fopen(fullfile(out, [fname '_acceptance.csv']), 'w');
    fprintf(fid, 'block,acceptance\n');
    acc = rec.record.AcceptanceRatio;
    for b = 1:numel(acc), fprintf(fid, '%d,%.6f\n', b, acc(b)); end
    fclose(fid);
    thin = str2double(getenv('SW_THIN'));
    if isnan(thin) || thin < 1, thin = 1; end
    fid = fopen(fullfile(out, [fname '_draws.csv']), 'w');
    fprintf(fid, 'block,iter,logpost,%s\n', strjoin(key', ','));
    nblck = options_.mh_nblck;
    for b = 1:nblck
        it = 0; k = 1;
        while true
            f = fullfile(mhdir, sprintf('%s_mh%d_blck%d.mat', fname, k, b));
            if ~exist(f, 'file'), break; end
            s = load(f, 'x2', 'logpo2');
            for r = 1:size(s.x2, 1)
                it = it + 1;
                if mod(it, thin) ~= 0, continue; end
                fprintf(fid, '%d,%d,%.17g', b, it, s.logpo2(r));
                fprintf(fid, ',%.17g', s.x2(r, :));
                fprintf(fid, '\n');
            end
            k = k + 1;
        end
    end
    fclose(fid);
end

% ---- Dynare's log-likelihood and log-prior at given points --------------
pts_file = fullfile(out, 'points.csv');
if exist(pts_file, 'file')
    % plain parser (works in MATLAB and Octave): header "point,<key>,..."
    fin = fopen(pts_file, 'r');
    hdr = strsplit(strtrim(fgetl(fin)), ',');
    rows = {};
    while true
        ln = fgetl(fin);
        if ~ischar(ln), break; end
        if ~isempty(strtrim(ln)), rows{end+1} = strsplit(strtrim(ln), ','); end
    end
    fclose(fin);
    col = zeros(nx, 1);
    for i = 1:nx
        j = find(strcmp(hdr, key{i}));
        if numel(j) ~= 1, error('points.csv has no column %s', key{i}); end
        col(i) = j;
    end
    fid = fopen(fullfile(out, [fname '_points_eval.csv']), 'w');
    fprintf(fid, 'point,loglik,logprior,logpost\n');
    for r = 1:numel(rows)
        x = zeros(nx, 1);
        for i = 1:nx, x(i) = str2double(rows{r}{col(i)}); end
        ll = evaluate_likelihood(x, M_, estim_params_, oo_, options_, bayestopt_);
        lp = priordens(x, bayestopt_.pshape, bayestopt_.p6, bayestopt_.p7, ...
                       bayestopt_.p3, bayestopt_.p4);
        fprintf(fid, '%s,%.12f,%.12f,%.12f\n', rows{r}{1}, ll, lp, ll + lp);
    end
    fclose(fid);
end
fprintf('Exported %s to %s\n', fname, out);
