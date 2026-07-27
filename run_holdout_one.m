here = fileparts(mfilename('fullpath'));
addpath(genpath(here));
P = repro_paths();
review_dir = P.outdir;
DATA_FILE  = P.data_file;

cfg = P.cfg;

valid_pids = P.cohort;

N = get_patient_count(DATA_FILE);

t0 = tic;
out = run_cohort_holdout(cfg, q_one, DATA_FILE, N);

models = {'pop','k1d','m2d','loglin3d','van4d'};
for m = models
    IN.(m{1})  = mean(out.mae_in.(m{1})(valid_pids),  'omitnan');
    OUT.(m{1}) = mean(out.mae_out.(m{1})(valid_pids), 'omitnan');
end

fprintf('q=%.0e  IN : 1D %.3f 2D %.3f 3D %.3f 4D %.3f pop %.3f\n', q_one, ...
    IN.k1d, IN.m2d, IN.loglin3d, IN.van4d, IN.pop);
fprintf('q=%.0e  OUT: 1D %.3f 2D %.3f 3D %.3f 4D %.3f pop %.3f   (%.0f s)\n', q_one, ...
    OUT.k1d, OUT.m2d, OUT.loglin3d, OUT.van4d, OUT.pop, toc(t0));

point = struct('q', q_one, 'IN', IN, 'OUT', OUT, ...
               'mae_in', out.mae_in, 'mae_out', out.mae_out, 'valid_pids', valid_pids);
save(fullfile(review_dir, sprintf('holdout_q%s.mat', tag)), 'point');
fprintf('saved holdout_q%s.mat\n', tag);
