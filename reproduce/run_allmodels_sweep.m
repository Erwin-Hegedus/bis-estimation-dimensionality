%RUN_ALLMODELS_SWEEP  In-sample MAE and parameter roughness against process noise.
%   Sweeps q for all five models and plots MAE against trajectory roughness, to
%   show whether lower in-sample error is bought by parameter thrashing.

addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
P = repro_paths();
cfg        = P.cfg;
DATA_FILE  = P.data_file;
valid_pids = P.cohort;
review_dir = P.outdir;
save_pids  = [1, 105, 150];

q_grid = [0, 1e-4, 1e-3, 1e-2, 1e-1, 1];
N  = get_patient_count(DATA_FILE);
nq = numel(q_grid);

models = {'k1d','m2d','loglin3d','van4d'};
res = cell(nq,1);
mae_med  = struct(); rough_med = struct();
for m = models, mae_med.(m{1}) = nan(1,nq); rough_med.(m{1}) = nan(1,nq); end
mae_pop = nan(1,nq);

for j = 1:nq
    fprintf('\n===== all-model sweep %d/%d: q = %.0e =====\n', j, nq, q_grid(j));
    t0 = tic;
    out = run_cohort_allmodels(cfg, q_grid(j), DATA_FILE, N, save_pids);
    res{j} = out;
    mae_pop(j) = mean(out.mae.pop(valid_pids), 'omitnan');
    for m = models
        mae_med.(m{1})(j)   = mean(out.mae.(m{1})(valid_pids), 'omitnan');
        rough_med.(m{1})(j) = median(out.rough.(m{1})(valid_pids), 'omitnan');
    end
    fprintf('q=%.0e  MAE 1D %.3f  2D %.3f  3D %.3f  4D %.3f  | rough 4D %.4f  in %.0f s\n', ...
        q_grid(j), mae_med.k1d(j), mae_med.m2d(j), mae_med.loglin3d(j), mae_med.van4d(j), ...
        rough_med.van4d(j), toc(t0));
    sweep = struct('q_grid', q_grid, 'valid_pids', valid_pids, 'mae_med', mae_med, ...
                   'rough_med', rough_med, 'mae_pop', mae_pop, 'res', {res});
    save(fullfile(review_dir, 'allmodels_sweep.mat'), 'sweep', '-v7.3');
end

fprintf('\n=== in-sample MAE vs q (cohort mean) ===\n');
for j = 1:nq
    fprintf('  q=%.0e  1D %.3f  2D %.3f  3D %.3f  4D %.3f  pop %.3f\n', q_grid(j), ...
        mae_med.k1d(j), mae_med.m2d(j), mae_med.loglin3d(j), mae_med.van4d(j), mae_pop(j));
end

x = max(q_grid, 1e-5);
col = struct('k1d',[0.2 0.4 0.7],'m2d',[0.2 0.6 0.3],'loglin3d',[0.7 0.5 0.1],'van4d',[0.8 0.3 0.2]);
lab = struct('k1d','1D','m2d','2D','loglin3d','3D','van4d','4D');

fig = figure('Color','w','Position',[60 60 1200 480]);
subplot(1,2,1); hold on;
for m = models
    semilogx(x, mae_med.(m{1}), '-o', 'LineWidth', 2, 'Color', col.(m{1}), 'DisplayName', lab.(m{1}));
end
set(gca,'XScale','log','FontSize',12); grid on;
xlabel('Process noise q  (\times initial covariance)', 'FontSize', 14);
ylabel('In-sample cohort MAE (BIS)', 'FontSize', 14);
title('(a) MAE vs adaptation, all models', 'FontSize', 15, 'FontWeight', 'bold');
legend('Location','best','FontSize',12);

subplot(1,2,2); hold on;
for m = models
    semilogx(x, rough_med.(m{1}), '-o', 'LineWidth', 2, 'Color', col.(m{1}), 'DisplayName', lab.(m{1}));
end
set(gca,'XScale','log','FontSize',12); grid on;
xlabel('Process noise q  (\times initial covariance)', 'FontSize', 14);
ylabel('Parameter trajectory roughness', 'FontSize', 14);
title('(b) Overfitting signature (roughness)', 'FontSize', 15, 'FontWeight', 'bold');
legend('Location','best','FontSize',12);

sgtitle('All-model Q sweep (N=209): does adaptation lower MAE by overfitting?', 'FontSize', 16, 'FontWeight', 'bold');
exportgraphics(fig, fullfile(review_dir, 'figure_allmodels_sweep.png'), 'Resolution', 300);
fprintf('Saved allmodels_sweep.mat and figure_allmodels_sweep.png\n');
