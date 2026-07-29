%RUN_FIM_RANK  Effective rank of the FIM across the cohort, at the estimator's
%   own threshold (ident_eigenvalue_ratio = 0.01). The filter estimates
%   log-parameters, so the spectrum is dimensionless and does not depend on the
%   units the parameters are carried in.
%
%   Both accumulations are reported:
%     FIM_cum    undiscounted sum over the whole case
%     FIM_final  the online discounted FIM at the last sample; with
%                lambda_F = 0.995 its memory is 1/(1-lambda_F) = 200 samples,
%                so it characterises only the final minutes of the case

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

sources = {'FIM_cum', 'FIM_final'};
out = struct();
missing = 0;

for s = 1:numel(sources)
    fld = sources{s};
    rank_raw = nan(nc,1); lam = nan(nc,4);
    share1 = nan(nc,1); share12 = nan(nc,1);
    gap12 = nan(nc,1); gap23 = nan(nc,1);

    for c = 1:nc
        i = coh(c);
        if ~isfield(results.raw(i), fld) || isempty(results.raw(i).(fld))
            if s == 1, missing = missing + 1; end
            continue;
        end
        v = sort(real(eig(results.raw(i).(fld))), 'descend');
        lam(c,:) = v.';
        rank_raw(c) = sum(v > RHO * v(1));
        share1(c)  = v(1) / sum(v);
        share12(c) = sum(v(1:2)) / sum(v);
        if v(1) > 0, gap12(c) = v(2) / v(1); end
        if v(2) > 0, gap23(c) = v(3) / v(2); end
    end

    out.(fld) = struct('rank_raw', rank_raw, 'lam', lam, 'share1', share1, ...
                       'share12', share12, 'gap12', gap12, 'gap23', gap23);

    fprintf('\n===== %s =====\n', fld);
    fprintf('EFFECTIVE RANK at the code''s own threshold (lambda_i > %.2g * lambda_1)\n', RHO);
    for k = 1:4
        n = sum(rank_raw == k);
        fprintf('  %d direction(s): %3d / %3d  (%.0f%%)\n', k, n, nc, 100*n/nc);
    end
    fprintf('SPECTRAL SHAPE (median over cohort)\n');
    fprintf('  lambda_1 share of total information : %.1f%%\n', 100*median(share1,'omitnan'));
    fprintf('  lambda_1+lambda_2 share             : %.2f%%\n', 100*median(share12,'omitnan'));
    fprintf('  lambda_2 / lambda_1                 : %.4g\n', median(gap12,'omitnan'));
    fprintf('  lambda_3 / lambda_2                 : %.4g\n', median(gap23,'omitnan'));
    fprintf('  median eigenvalues                  : [%.4g %.4g %.4g %.4g]\n', median(lam,1,'omitnan'));
end

fprintf('\ncohort %d, patients missing FIM_cum: %d\n', nc, missing);

prov = provenance_stamp(cfg, P.data_file);
save(fullfile(P.outdir,'fim_rank.mat'), 'out', 'coh', 'RHO', 'prov');
fprintf('\nSaved %s\n', fullfile(P.outdir,'fim_rank.mat'));
