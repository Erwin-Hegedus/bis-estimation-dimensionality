%REGENERATE_PAPER_FIGURES  Rebuild every manuscript figure from the cached run.
%
%   The cohort loop in run_analysis.m takes hours and its output is already on
%   disk, so styling work reads results/bis_analysis_results_v6_0.mat instead of
%   re-running it. Nothing here recomputes an estimate; if a number changes,
%   something is wrong.
%
%   Figure names on the left are what the .tex expects; the copy step at the end
%   pushes them into 'paper jbhi revised/' under those names.

src = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(src));
cd(src);

fprintf('Loading cached results...\n');
S   = load(fullfile(src, 'results', 'bis_analysis_results_v6_0.mat'), ...
           'results', 'cfg', 'eq_info', 'pid_test');
results = S.results;
cfg     = S.cfg;
eq_info = S.eq_info;
clear S;

fig_dir = fullfile(src, 'results');
close all;

%% ---- figures the manuscript includes ---------------------------------
D = load(fullfile(src, 'patientDataFinal_auto.mat'), 'patientDataFinal');
generate_figure1_protocol(D.patientDataFinal, fig_dir);
clear D;

generate_figure2_mae_boxplots(results, fig_dir);
generate_figure4_case_study(results, 107, cfg, fig_dir);
generate_figure5_monte_carlo(eq_info, cfg, fig_dir);
generate_combined_identifiability_figure(results, cfg, fig_dir, 105);

plot_parameter_clustering_boxplots(results, fig_dir);
plot_population_stats_boxplot(results, fig_dir);

vec   = build_vec_from_results(results);
cases = pick_case_studies(vec, results);
plot_case_studies(cases, results, cfg, fig_dir);
plot_case_study_params(cases, results, cfg, fig_dir);

%% ---- supporting figures kept in results/ ------------------------------
generate_figure3_dimensionality(results, fig_dir);
generate_figure6_2d_drift(results, fig_dir);
plot_1d_vs_2d(results, 105, fig_dir);

close all;

%% ---- push into the manuscript directory -------------------------------
paper_dir = fullfile(fileparts(src), 'paper jbhi revised');
copy_map = { ...
    'figure1_protocol.png',                  'figure1_protocol.png'; ...
    'figure2_mae_boxplots.png',              'figure2_mae_boxplots.png'; ...
    'figure4_case_study.png',                'figure4_case_study_1d2d.png'; ...
    'figure5_monte_carlo.png',               'figure5_monte_carlo.png'; ...
    'figure7_combined_identifiability.png',  'figure7_combined_ident.png'; ...
    'figure_4d_params.png',                  'figure_4d_params.png'; ...
    'figure_cases_comparison.png',           'figure_cases_comparison2.png'; ...
    'figure_param_boxplots_1d2d.png',        'figure_param_boxplots_1d2d.png'; ...
    'figure_population_boxplot.png',         'figure_population_boxplot.png'};

if exist(paper_dir, 'dir')
    for k = 1:size(copy_map, 1)
        copyfile(fullfile(fig_dir, copy_map{k,1}), fullfile(paper_dir, copy_map{k,2}));
        fprintf('  -> %s\n', copy_map{k,2});
    end
else
    warning('regen:paperdir', 'Manuscript directory not found: %s', paper_dir);
end

fprintf('\nDone.\n');
