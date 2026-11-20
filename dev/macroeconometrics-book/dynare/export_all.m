% Export several Dynare runs with export_dynare.m (MATLAB or Octave).
% SW_EXPORT_LIST: file with one run per line: <run_dir>|<fname>|<results.mat>|<out_dir>|<thin>
% SW_SCRIPT_DIR:  this folder (used when mfilename cannot tell)
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = getenv('SW_SCRIPT_DIR'); end
list_file = getenv('SW_EXPORT_LIST');
fin = fopen(list_file, 'r');
runs = {};
while true
    ln = fgetl(fin);
    if ~ischar(ln), break; end
    if ~isempty(strtrim(ln)), runs{end+1} = strsplit(strtrim(ln), '|'); end
end
fclose(fin);
for r_ = 1:numel(runs)
    spec_ = runs{r_};
    setenv('SW_RUN_DIR', spec_{1}); setenv('SW_FNAME', spec_{2});
    setenv('SW_RESULTS', spec_{3}); setenv('SW_OUT', spec_{4}); setenv('SW_THIN', spec_{5});
    clear evaluate_likelihood dsge_likelihood
    run(fullfile(script_dir, 'export_dynare.m'));
    clearvars -except runs r_ script_dir
end
