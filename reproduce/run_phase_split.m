%RUN_PHASE_SPLIT  Induction vs maintenance MAE decomposition, computed offline
%   from results.raw. First 10 samples skipped, matching run_analysis.
%     DEF A  induction = start -> first sample with BIS < 60
%     DEF B  induction = start -> start + 10 min

addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
P = repro_paths();

canonical = fullfile(P.repo, 'results', 'bis_analysis_results_v6_0.mat');
if ~exist(canonical, 'file')
    error('repro:missingCanonical', '%s not found; run run_analysis.m first.', canonical);
end
fprintf('loading %s ...\n', canonical);
S = load(canonical, 'results');
results = S.results; clear S;

models = {'pop','k1d','m2d','loglin3d','van4d','gre4d'};
field  = {'pred_pop','pred_kscale','pred_2d','pred_loglin','pred_van','pred_gre'};
label  = {'population','1D','2D','3D','4D Bouillon','4D Greco'};

N0 = 10;              % initialization samples skipped
IND_SECONDS = 600;    % DEF B induction window

nm = numel(models);
cohort = P.cohort;
nc = numel(cohort);
mae = struct();
for d = {'A','B'}
    for m = models
        mae.(d{1}).ind.(m{1})  = nan(nc,1);
        mae.(d{1}).main.(m{1}) = nan(nc,1);
    end
end
mae.whole = struct();
for m = models, mae.whole.(m{1}) = nan(nc,1); end
ind_share = nan(nc,2);
caselen = nan(nc,1);

for c = 1:nc
    i = cohort(c);
    r = results.raw(i);
    if isempty(r.bis) || isempty(r.time), continue; end
    bis = r.bis(:); t = r.time(:);
    n = numel(bis);
    caselen(c) = n;
    keep = (N0:n)';

    % boundary indices
    below = find(bis < 60);
    below = below(below >= N0);
    iA = []; if ~isempty(below), iA = below(1); end
    iB = find(t <= t(1) + IND_SECONDS, 1, 'last');

    for k = 1:nm
        p = r.(field{k});
        if isempty(p), continue; end
        p = p(:);
        mae.whole.(models{k})(c) = mean(abs(p(keep) - bis(keep)), 'omitnan');

        for d = {'A','B'}
            if strcmp(d{1},'A'), ib = iA; else, ib = iB; end
            if isempty(ib) || ib <= N0 || ib >= n, continue; end
            vi = (N0:ib)';
            vm = ((ib+1):n)';
            mae.(d{1}).ind.(models{k})(c)  = mean(abs(p(vi) - bis(vi)), 'omitnan');
            mae.(d{1}).main.(models{k})(c) = mean(abs(p(vm) - bis(vm)), 'omitnan');
        end
    end
    if ~isempty(iA), ind_share(c,1) = iA / n; end
    if ~isempty(iB), ind_share(c,2) = iB / n; end
end

fprintf('\nmedian case length %.0f samples (%.0f min)\n', ...
    median(caselen,'omitnan'), median(caselen,'omitnan')/60);
fprintf('induction share: DEF A median %.1f%%   DEF B median %.1f%%\n', ...
    100*median(ind_share(:,1),'omitnan'), 100*median(ind_share(:,2),'omitnan'));

for d = {'B','A'}
    usable = sum(~isnan(mae.(d{1}).ind.van4d));
    fprintf('\n=== DEF %s  (%d/%d patients usable) ===\n', d{1}, usable, nc);
    fprintf('%-14s %10s %12s %11s\n', 'model', 'INDUCTION', 'MAINTENANCE', 'whole case');
    for k = 1:nm
        fprintf('%-14s %10.2f %12.2f %11.2f\n', label{k}, ...
            mean(mae.(d{1}).ind.(models{k}), 'omitnan'), ...
            mean(mae.(d{1}).main.(models{k}), 'omitnan'), ...
            mean(mae.whole.(models{k}), 'omitnan'));
    end
    adaptive = models(~strcmp(models,'pop'));
    v = nan(numel(adaptive),1);
    for k = 1:numel(adaptive), v(k) = mean(mae.(d{1}).main.(adaptive{k}), 'omitnan'); end
    fprintf('maintenance spread, adaptive models: %.2f BIS\n', max(v) - min(v));
    v3 = v(~strcmp(adaptive','loglin3d'));
    fprintf('maintenance spread, excluding 3D:    %.2f BIS\n', max(v3) - min(v3));
end

prov = provenance_stamp(P.cfg, P.data_file);
save(fullfile(P.outdir, 'phase_split.mat'), 'mae', 'ind_share', 'caselen', 'cohort', 'models', 'label', 'prov');
fprintf('\nSaved %s\n', fullfile(P.outdir, 'phase_split.mat'));
