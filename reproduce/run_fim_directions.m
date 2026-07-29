%RUN_FIM_DIRECTIONS  Which parameter directions the BIS data identifies.
%   The filter estimates log-parameters, so the accumulated FIM is already
%   dimensionless and no rescaling is applied. The alternative scalings this
%   script used to report cannot be recovered from it: the accumulation runs at
%   the theta(t) of each sample, so no single W inverts it.
%
%   Set MODEL to 'van' (default) or 'gre' before running.

addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
P = repro_paths();
if ~exist('MODEL','var'), MODEL = 'van'; end

canonical = fullfile(P.repo, 'results', 'bis_analysis_results_v6_0.mat');
if ~exist(canonical, 'file')
    error('repro:missingCanonical', '%s not found; run run_analysis.m first.', canonical);
end
S = load(canonical, 'results', 'cfg');
results = S.results; cfg = S.cfg; clear S;

coh = P.cohort;
switch MODEL
    case 'van', fimField = 'FIM_final';     xField = 'Xhist_van';
    case 'gre', fimField = 'FIM_final_gre'; xField = 'Xhist_gre';
    otherwise,  error('MODEL must be ''van'' or ''gre''');
end
if ~isfield(results.raw, fimField)
    error('repro:missingFIM', '%s not in results; rerun run_analysis.m.', fimField);
end
nm = {'C50P','C50R','gamma','synergy'};

V1 = []; V2 = []; L = []; nb = 0;
for c = 1:numel(coh)
    i = coh(c);
    F = results.raw(i).(fimField);
    if isempty(F) || any(~isfinite(F(:))), continue; end

    X = results.raw(i).(xField);
    k = find(any(~isnan(X),1), 1, 'last');
    if isempty(k), continue; end
    th = X(:,k);
    if any(abs(th-cfg.lb(:))<1e-9 | abs(th-cfg.ub(:))<1e-9), nb = nb + 1; end

    A = (F+F')/2;
    [V,D] = eig(A); d = real(diag(D));
    [d,o] = sort(d,'descend'); V = real(V(:,o));
    al = @(v) v * sign(v(find(abs(v)==max(abs(v)),1)));
    V1(:,end+1) = al(V(:,1)); V2(:,end+1) = al(V(:,2)); L(:,end+1) = d;
end

n = size(V1,2);
[~,d1] = max(abs(V1),[],1); [~,d2] = max(abs(V2),[],1);
rk = sum(L > 0.01 * L(1,:), 1);

fprintf('\n===== [%s] log-FIM, fractional parameter directions  (n=%d) =====\n', MODEL, n);
fprintf('  fitted parameter on a bound: %d/%d patients\n', nb, n);
fprintf('  direction 1 led by: '); for k=1:4, fprintf('%s %.0f%%  ', nm{k}, 100*mean(d1==k)); end; fprintf('\n');
fprintf('  direction 2 led by: '); for k=1:4, fprintf('%s %.0f%%  ', nm{k}, 100*mean(d2==k)); end; fprintf('\n');
fprintf('  median dir1 = [%+.2f %+.2f %+.2f %+.2f]\n', median(V1,2));
fprintf('  median dir2 = [%+.2f %+.2f %+.2f %+.2f]\n', median(V2,2));
dd = abs(V1(1,:)) - abs(V1(2,:));
fprintf('  |C50P|-|C50R| in dir1: median %+.3f  IQR [%+.3f %+.3f]\n', ...
    median(dd), prctile(dd,25), prctile(dd,75));
fprintf('  effective rank: 1 in %.0f%%, 2 in %.0f%%, 3 in %.0f%%, 4 in %.0f%%\n', ...
    100*mean(rk==1), 100*mean(rk==2), 100*mean(rk==3), 100*mean(rk==4));

out = struct('V1',V1,'V2',V2,'lambda',L,'rank',rk,'nb',nb, ...
             'label','log-FIM, fractional parameter directions');

prov = provenance_stamp(cfg, P.data_file);
save(fullfile(P.outdir,['fim_directions_' MODEL '.mat']), 'out', 'coh', 'nm', 'MODEL', 'prov');
fprintf('\nSaved %s\n', fullfile(P.outdir,['fim_directions_' MODEL '.mat']));
