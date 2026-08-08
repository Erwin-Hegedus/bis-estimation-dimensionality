src = '/Users/erwinhegedus/Library/CloudStorage/Dropbox/Conferences/2025/JBHI/Source code';
addpath(genpath(src));
S = load(fullfile(src,'results','bis_analysis_results_v6_0.mat'),'results','cfg');
r = S.results; cfg = S.cfg;
coh = readmatrix(fullfile(src,'reproduce','cohort_caseids.csv')); coh = coh(:,1);

%% R3-2 : are CeP and CeR collinear over the cohort?
rho = nan(numel(coh),1); rho_log = nan(numel(coh),1); ratio_cv = nan(numel(coh),1);
for c = 1:numel(coh)
    a = r.raw(coh(c)).CeP_trajectory; b = r.raw(coh(c)).CeR_trajectory;
    if isempty(a), continue; end
    k0 = cfg.online.initialization_samples;
    a = a(k0:end); b = b(k0:end);
    ok = isfinite(a) & isfinite(b) & a > 0 & b > 0;
    if sum(ok) < 100, continue; end
    a = a(ok); b = b(ok);
    rho(c)     = corr(a, b);
    rho_log(c) = corr(log(a), log(b));
    q = (b*1000) ./ a;                 % remi:prop effect-site ratio
    ratio_cv(c) = std(q) / mean(q);
end
f = @(v) sprintf('median %.3f   IQR [%.3f, %.3f]   p10 %.3f   p90 %.3f', ...
    median(v,'omitnan'), prctile(v,25), prctile(v,75), prctile(v,10), prctile(v,90));
fprintf('=== R3-2  CeP vs CeR collinearity (N=%d) ===\n', sum(isfinite(rho)));
fprintf('  Pearson r          : %s\n', f(rho));
fprintf('  Pearson r on logs  : %s\n', f(rho_log));
fprintf('  frac with r > 0.9  : %.1f%%\n', 100*mean(rho > 0.9, 'omitnan'));
fprintf('  frac with r > 0.95 : %.1f%%\n', 100*mean(rho > 0.95,'omitnan'));
fprintf('  CV of CeR/CeP ratio: %s\n', f(ratio_cv));

%% R3-5 : is the leading FIM direction the shared-gain direction?
D = load(fullfile(src,'reproduce','output','fim_directions_van.mat'));
V1 = D.out.V1;                       % 4 x n, unit eigenvectors (log coords)
u  = [1;1;0;0]/sqrt(2);              % equal fractional change in C50P and C50R
al = abs(u' * V1);                   % |cos| with the shared-gain direction
c50 = sqrt(sum(V1(1:2,:).^2, 1));    % weight of dir1 on the two C50s
same_sign = sign(V1(1,:)) == sign(V1(2,:));
fprintf('\n=== R3-5  leading FIM direction vs the shared-gain direction (n=%d) ===\n', size(V1,2));
fprintf('  |cos| with [1,1,0,0]/sqrt2 : %s\n', f(al'));
fprintf('  angle to it (deg)          : %s\n', f(acosd(min(al,1))'));
fprintf('  C50P and C50R same sign    : %.1f%% of patients\n', 100*mean(same_sign));
fprintf('  weight of dir1 on the C50s : %s\n', f(c50'));
fprintf('  median dir1 (log coords)   : [%+.3f %+.3f %+.3f %+.3f]\n', median(V1,2));
