%RUN_HOLDOUT_SWEEP  Held-out evaluation of all models against process noise.
%   Adapts on the first half of each case, freezes, predicts the second half,
%   and reports in-sample and held-out cohort MAE at each q.

addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
P = repro_paths();
cfg        = P.cfg;
DATA_FILE  = P.data_file;
valid_pids = P.cohort;
review_dir = P.outdir;

% Set q_grid and out_tag before running to override these defaults.
if ~exist('q_grid', 'var'),  q_grid  = [0, 1e-4, 1e-3, 1e-2]; end
if ~exist('out_tag', 'var'), out_tag = ''; end
N  = get_patient_count(DATA_FILE);
nq = numel(q_grid);

models = {'pop','k1d','m2d','loglin3d','van4d','const'};
res = cell(nq,1);
for m = models
    IN.(m{1})  = nan(1,nq);
    OUT.(m{1}) = nan(1,nq);
end

for j = 1:nq
    fprintf('\n===== holdout sweep %d/%d: q = %.0e =====\n', j, nq, q_grid(j));
    t0 = tic;
    out = run_cohort_holdout(cfg, q_grid(j), DATA_FILE, N, valid_pids);
    res{j} = out;
    for m = models
        IN.(m{1})(j)  = mean(out.mae_in.(m{1})(valid_pids),  'omitnan');
        OUT.(m{1})(j) = mean(out.mae_out.(m{1})(valid_pids), 'omitnan');
    end
    fprintf('q=%.0e  IN : 1D %.3f 2D %.3f 3D %.3f 4D %.3f pop %.3f\n', q_grid(j), ...
        IN.k1d(j), IN.m2d(j), IN.loglin3d(j), IN.van4d(j), IN.pop(j));
    fprintf('q=%.0e  OUT: 1D %.3f 2D %.3f 3D %.3f 4D %.3f pop %.3f const %.3f  (%.0f s)\n', q_grid(j), ...
        OUT.k1d(j), OUT.m2d(j), OUT.loglin3d(j), OUT.van4d(j), OUT.pop(j), OUT.const(j), toc(t0));
    holdout = struct('q_grid', q_grid, 'valid_pids', valid_pids, 'IN', IN, 'OUT', OUT, 'res', {res});
    prov = provenance_stamp(cfg, DATA_FILE);
    save(fullfile(review_dir, ['holdout_sweep' out_tag '.mat']), 'holdout', 'prov', '-v7.3');
end

fprintf('\n=== SUMMARY: in-sample (first half, adapting) vs holdout (second half, frozen) ===\n');
for j = 1:nq
    fprintf('q=%.0e\n', q_grid(j));
    for m = {'k1d','m2d','loglin3d','van4d'}
        fprintf('   %-9s in %.3f  out %.3f  (gap %+.3f)\n', m{1}, IN.(m{1})(j), OUT.(m{1})(j), OUT.(m{1})(j)-IN.(m{1})(j));
    end
end

x = max(q_grid, 1e-5);
col = struct('k1d',[0.2 0.4 0.7],'m2d',[0.2 0.6 0.3],'loglin3d',[0.7 0.5 0.1],'van4d',[0.8 0.3 0.2]);
lab = struct('k1d','1D','m2d','2D','loglin3d','3D','van4d','4D');
mm = {'k1d','m2d','loglin3d','van4d'};

fig = figure('Color','w','Position',[60 60 1200 480]);
subplot(1,2,1); hold on;
for m = mm
    plot(x, IN.(m{1}), '-o', 'LineWidth', 2, 'Color', col.(m{1}), 'DisplayName', lab.(m{1}));
end
set(gca,'XScale','log','FontSize',12); grid on;
xlabel('Process noise q', 'FontSize', 14);
ylabel('MAE (BIS)', 'FontSize', 14);
title('(a) In-sample (adapting)', 'FontSize', 15, 'FontWeight', 'bold');
legend('Location','best','FontSize',12);

subplot(1,2,2); hold on;
for m = mm
    plot(x, OUT.(m{1}), '-o', 'LineWidth', 2, 'Color', col.(m{1}), 'DisplayName', lab.(m{1}));
end
plot(x, OUT.pop, 'k--', 'LineWidth', 2, 'DisplayName', 'population');
plot(x, OUT.const, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'DisplayName', 'constant (no drug data)');
set(gca,'XScale','log','FontSize',12); grid on;
xlabel('Process noise q', 'FontSize', 14);
ylabel('MAE (BIS)', 'FontSize', 14);
title('(b) Held-out second half (frozen)', 'FontSize', 15, 'FontWeight', 'bold');
legend('Location','best','FontSize',12);

sgtitle('Out-of-sample test: does lower in-sample error transfer?', 'FontSize', 16, 'FontWeight', 'bold');
export_figure_png(fig, fullfile(review_dir, ['figure_holdout_sweep' out_tag '.png']));
fprintf('Saved holdout_sweep.mat and figure_holdout_sweep.png\n');
