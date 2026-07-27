%RUN_FIM_RANK  Effective rank of the cumulative FIM across the cohort, at the
%   estimator's own threshold (ident_eigenvalue_ratio = 0.01). Spectrum is in
%   raw parameter coordinates, so it is unit-dependent.

addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
P = repro_paths();

canonical = fullfile(P.repo, 'results', 'bis_analysis_results_v6_0.mat');
if ~exist(canonical, 'file')
    error('repro:missingCanonical', '%s not found; run run_analysis.m first.', canonical);
end
fprintf('loading %s ...\n', canonical);
S = load(canonical, 'results', 'cfg');
results = S.results; cfg = S.cfg; clear S;

RHO = 0.01;                       % ekf/init_ekf_4d.m ident_eigenvalue_ratio
coh = P.cohort; nc = numel(coh);

rank_raw = nan(nc,1);
lam = nan(nc,4);
share1 = nan(nc,1);
share12 = nan(nc,1);
gap23 = nan(nc,1);
missing = 0;

for c = 1:nc
    i = coh(c);
    if ~isfield(results.raw(i), 'FIM_eig_hist') || isempty(results.raw(i).FIM_eig_hist)
        missing = missing + 1; continue;
    end
    E = results.raw(i).FIM_eig_hist;
    k = find(any(~isnan(E),1), 1, 'last');
    if isempty(k), continue; end
    v = sort(real(E(:,k)), 'descend');
    lam(c,:) = v.';
    rank_raw(c) = sum(v > RHO * v(1));
    share1(c)  = v(1) / sum(v);
    share12(c) = sum(v(1:2)) / sum(v);
    if v(3) > 0, gap23(c) = v(2) / v(3); end
end

fprintf('\ncohort %d, patients missing FIM_eig_hist: %d\n', nc, missing);
fprintf('\nEFFECTIVE RANK at the code''s own threshold (lambda_i > %.2g * lambda_1)\n', RHO);
for k = 1:4
    n = sum(rank_raw == k);
    fprintf('  %d direction(s): %3d / %3d  (%.0f%%)\n', k, n, nc, 100*n/nc);
end
fprintf('\nSPECTRAL SHAPE (median over cohort)\n');
fprintf('  lambda_1 share of total information : %.1f%%\n', 100*median(share1,'omitnan'));
fprintf('  lambda_1+lambda_2 share             : %.2f%%\n', 100*median(share12,'omitnan'));
fprintf('  lambda_2 / lambda_3                 : %.0f\n', median(gap23,'omitnan'));
fprintf('  median eigenvalues                  : [%.4g %.4g %.4g %.4g]\n', median(lam,1,'omitnan'));

save(fullfile(P.outdir,'fim_rank.mat'), 'rank_raw', 'lam', 'share1', 'share12', 'gap23', 'coh', 'RHO');
fprintf('\nSaved %s\n', fullfile(P.outdir,'fim_rank.mat'));
