%% Path Setup
addpath(fullfile(fileparts(mfilename('fullpath')), 'core'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'ekf'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'physiology'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'processor'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'data'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'analysis'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'figures'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'diagnostic_plots'));

clearvars -except -regexp ^OCTAVE_VERSION$; clc; close all;
cfg = struct();
DATA_FILE = 'patientDataFinal_auto.mat';
if ~exist(DATA_FILE, 'file'); error('Data file %s not found.', DATA_FILE); end

% ================= CORE CONFIG =====================
cfg.ke0P = 0.456;  % Schnider population (DO NOT inflate)
cfg.ke0R = 0.595;  % Minto population
cfg.BISmin_fixed = 30;
cfg.E0_fixed = 93;
cfg.bis_delay = 25;  % samples at 1 Hz, BIS monitor smoothing lag

% Gate settings
cfg.gate = struct('window',60,'forget',0.98,'temperature',2.0,'min_weight',0.00);

cfg.pop_virtual_sigma = [0.50, 0.60, 0.30, 0.30];
cfg.decision_point = struct( ...
    'evaluation_start',300,'evaluation_duration',300,'min_samples',100, ...
    'improvement_threshold',1.0,'significance_level',0.05, ...
    'effect_size_threshold',0.3,'consistency_threshold',3.0);

cfg.R_population = diag([0.1^2, 0.15^2, 0.15^2, 0.08^2]);
cfg.population_params_van = [3.5, 5.0, 2.8, 0.8];
cfg.population_params_gre = [3.2, 4.5, 2.8, 0.9];

% Physiological Bounds
cfg.lb = [1.0,  2.0,  1.0,  0.3]; 
cfg.ub = [12.0, 40.0, 6.0, 1.8];

cfg.P_min = [0.01, 0.02, 0.005, 0.002];
cfg.P_max = [4.0, 9.0, 1.5, 0.5];

cfg.param_rate_max = [0.02, 0.05, 0.01, 0.005];

cfg.fim_condition_threshold = 100;
cfg.fim_eigenvalue_floor = 1e-6;
cfg.fim_forgetting = 0.995;

cfg.R_base = 400;
cfg.R_disequilibrium_factor = 0;

cfg.pk.prop = struct('V1',4.27,'V2',18.9,'V3',238,'CL',1.89,'Q2',1.29,'Q3',0.836,'input_conc',20);
cfg.pk.remi = struct('V1',5.1,'V2',9.82,'V3',5.42,'CL',2.6,'Q2',2.05,'Q3',0.076,'input_conc',0.02);

cfg.target_band = [40,60];
cfg.min_Ce_for_learning = 0.1;
cfg.online = struct('initialization_samples',10);
cfg.artifact = struct('bis_min',0,'bis_max',100,'bis_jump_thresh',25,'bis_rate_thresh',15);
% Fisher Information Thresholds
cfg.fisher_threshold = [5; 5; 4; 3];
cfg.identifiability_threshold = 0.15;
cfg.R_floor = 80;

% ke0 Estimator Config
cfg.ke0_estimator = struct('lambda', 0.1, 'Q_ke0', (0.002/60)^2, 'adapt_window_min', 15);

fprintf('=== ONLINE BIS ANALYSIS V6.0 (Convergent Estimator) ===\n');
fprintf('Key Changes: Q=0, FIM Identifiability, Covariance Bounds, Rate Limiting\n');
fprintf('Models: 0D (Pop), 1D (k_scale), 2D (kP,kR), 3D (LogLin), 4D (Vanluchene)\n');

%% ======================= Main Processing Loop ===========================
results_dir = 'results';
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

N = get_patient_count(DATA_FILE);
results = init_results_struct(N);

cohort = readmatrix(fullfile('reproduce', 'cohort_caseids.csv'));
cohort = cohort(:, 1);
fprintf('Cohort: %d cases from reproduce/cohort_caseids.csv\n', numel(cohort));

for c = 1:numel(cohort)
    i = cohort(c);
    pid = i;
    fprintf('\nPatient %3d (%2d/%2d): ', pid, c, numel(cohort));
    try
        [bis, remi_rate, prop_rate, time] = load_patient_record(pid, DATA_FILE);
        n = min([numel(bis), numel(remi_rate), numel(prop_rate), numel(time)]);

        if n < 900
            fprintf('SKIP: Short case (<15 mins)\n'); continue;
        end

        bis = bis(1:n);
        remi_rate_raw = remi_rate(1:n);
        prop_rate_raw = prop_rate(1:n);
        time = time(1:n);


        demo = get_patient_demo_from_mat(pid, DATA_FILE);
        [pk_prop, pk_remi] = pk_from_demographics_schnider_minto(demo, cfg);
        [prop_rate, remi_rate, bis, time, shift_val, n, sync_msg] = ...
            align_bis_to_infusion( ...
            bis, prop_rate_raw, remi_rate_raw, time, 0.3 / 60 );

        fprintf('  -> Ce bootstrap alignment: %s (Shift %.1fs, n=%d)\n', ...
            sync_msg, shift_val, n);
        processor = init_processor(cfg, time(1), pk_prop, pk_remi, demo);

        pred_van = nan(n,1); pred_gre = nan(n,1); pred_pop = nan(n,1);
        pred_kscale = nan(n,1);
        pred_loglin = nan(n,1);
        Xhist_van = nan(4,n); Xhist_gre = nan(4,n);
        Xhist_loglin = nan(3,n);
        RSE_van = nan(4,n); RSE_gre = nan(4,n);
        Phist_van = nan(4,n);
        FIM_cond_hist = nan(n,1);
        FIM_eig_hist = nan(4,n);
        CeP_trajectory = nan(n,1); CeR_trajectory = nan(n,1);
        ke0_trajectory = nan(n,1);
        tau_d_trajectory = nan(n,1);
        BISmin_trajectory = nan(n,1);
        E0_trajectory = nan(n,1);
        k_trajectory = nan(n,1);
        pred_2d_fim = nan(n,1);
        kP_fim_trajectory = nan(n,1);
        kR_fim_trajectory = nan(n,1);
        kP_var_trajectory = nan(n,1);
        kR_var_trajectory = nan(n,1);

        
        for k = 1:n
            [pred_van(k), pred_gre(k), pred_kscale(k), pred_loglin(k), pred_2d_fim(k), processor] = process_sample( ...
                processor, time(k), bis(k), prop_rate(k), remi_rate(k));
            
            if processor.initialized
                CeP = processor.effect_site_P.Ce_delayed_output;
                CeR = processor.effect_site_R.Ce_delayed_output;
                
                CeP_trajectory(k) = CeP;
                CeR_trajectory(k) = CeR;
                ke0_trajectory(k) = processor.effect_site_P.ke0 * 60;
                tau_d_trajectory(k) = processor.effect_site_P.tau_d;
                
                pred_pop(k) = predict_bis_population(cfg.population_params_van, ...
                                                    CeP, CeR, cfg.E0_fixed, ...
                                                    'vanluchene', cfg.BISmin_fixed);

                k_trajectory(k) = processor.ekf_k.k;
                kP_fim_trajectory(k) = processor.ekf_2d_fim.current_params(1);
                kR_fim_trajectory(k) = processor.ekf_2d_fim.current_params(2);
                kP_var_trajectory(k) = processor.ekf_2d_fim.P(1,1);
                kR_var_trajectory(k) = processor.ekf_2d_fim.P(2,2);
                Xhist_loglin(:,k) = processor.ekf_loglin.current_params;
            else
                CeP_trajectory(k) = 0; CeR_trajectory(k) = 0;
                ke0_trajectory(k) = NaN; tau_d_trajectory(k) = NaN;
                pred_pop(k) = NaN;
                k_trajectory(k) = 1.0;
            end
            
            P_diag = diag(processor.ekf_van.P);
            if numel(P_diag) == 3
                Phist_van(:,k) = [P_diag; 0];
            else
                Phist_van(:,k) = P_diag;
            end
            
            if isfield(processor.ekf_van, 'FIM_condition')
                FIM_cond_hist(k) = processor.ekf_van.FIM_condition;
            end
            if isfield(processor.ekf_van, 'FIM_eigenvalues')
                FIM_eig_hist(:,k) = processor.ekf_van.FIM_eigenvalues(:);
            end
            
            if ~isempty(processor.ekf_van.current_params)
                Xhist_van(:,k) = processor.ekf_van.current_params(:);
            end
            BISmin_trajectory(k) = processor.BISmin.x;
            E0_trajectory(k) = processor.E0.x;
            if isfield(processor.ekf_van,'RSE_current') && ~isempty(processor.ekf_van.RSE_current)
                RSE_van(:,k) = processor.ekf_van.RSE_current(:);
            end
            if ~isempty(processor.ekf_gre.current_params)
                Xhist_gre(:,k) = processor.ekf_gre.current_params(:);
            end
            if isfield(processor.ekf_gre,'RSE_current') && ~isempty(processor.ekf_gre.RSE_current)
                RSE_gre(:,k) = processor.ekf_gre.RSE_current(:);
            end
        end
        
        valid_idx = max(1, processor.online.initialization_samples):n;
        if numel(valid_idx) > 10
            mae_v = mean(abs(pred_van(valid_idx) - bis(valid_idx)), 'omitnan');
            mae_g = mean(abs(pred_gre(valid_idx) - bis(valid_idx)), 'omitnan');
            mae_p = mean(abs(pred_pop(valid_idx) - bis(valid_idx)), 'omitnan');
            mae_k = mean(abs(pred_kscale(valid_idx) - bis(valid_idx)), 'omitnan');
            mae_ll = mean(abs(pred_loglin(valid_idx) - bis(valid_idx)), 'omitnan');
            mae_2d_fim = mean(abs(pred_2d_fim(valid_idx) - bis(valid_idx)), 'omitnan'); 
        else
            mae_v = NaN; mae_g = NaN; mae_p = NaN; mae_k = NaN; mae_ll = NaN; mae_2d_fim = NaN;
        end
        
        fprintf('OK (Pop=%.2f | K1D=%.2f | 2D=%.2f | LL3D=%.2f | Van4D=%.2f | Gre4D=%.2f)\n', ...
            mae_p, mae_k, mae_2d_fim, mae_ll, mae_v, mae_g);
        
        % Store demographics for Table 1
        results.demographics(i).Age = demo.Age;
        results.demographics(i).Weight = demo.Wt_kg;
        results.demographics(i).Height = demo.Ht_cm;
        results.demographics(i).Sex = demo.Sex;
        results.demographics(i).Duration = n / 60;  % minutes
        
        results.patient_id(i) = pid;
        results.raw(i).time = time(:);
        results.raw(i).bis = bis(:);
        results.raw(i).pred_van = pred_van(:);
        results.raw(i).pred_gre = pred_gre(:);
        results.raw(i).pred_pop = pred_pop(:);
        results.raw(i).pred_kscale = pred_kscale(:);
        results.raw(i).pred_loglin = pred_loglin(:);
        results.raw(i).pred_2d = pred_2d_fim(:); 
        results.raw(i).Xhist_van = Xhist_van;
        results.raw(i).Xhist_gre = Xhist_gre;
        results.raw(i).Xhist_loglin = Xhist_loglin;
        results.raw(i).RSE_van = RSE_van;
        results.raw(i).RSE_gre = RSE_gre;
        results.raw(i).CeP_trajectory = CeP_trajectory(:);
        results.raw(i).CeR_trajectory = CeR_trajectory(:);
        results.raw(i).ke0_trajectory = ke0_trajectory(:);
        results.raw(i).tau_d_trajectory = tau_d_trajectory(:);
        results.raw(i).Phist_van = Phist_van;
        results.raw(i).FIM_cond_hist = FIM_cond_hist(:);
        results.raw(i).FIM_eig_hist = FIM_eig_hist;
        results.raw(i).FIM_final = processor.ekf_van.FIM;
        results.raw(i).FIM_final_gre = processor.ekf_gre.FIM;
        results.raw(i).BISmin_trajectory = BISmin_trajectory(:);
        results.raw(i).E0_trajectory = E0_trajectory(:);
        results.raw(i).k_trajectory = k_trajectory(:);
        results.raw(i).kP_trajectory = kP_fim_trajectory(:);
        results.raw(i).kR_trajectory = kR_fim_trajectory(:);
        results.raw(i).pred_2d_fim = pred_2d_fim(:);
        results.raw(i).kP_fim_trajectory = kP_fim_trajectory(:);
        results.raw(i).kR_fim_trajectory = kR_fim_trajectory(:);
        results.raw(i).kP_hist = kP_fim_trajectory(:);
        results.raw(i).kR_hist = kR_fim_trajectory(:);
        results.raw(i).kP_var  = kP_var_trajectory(:);
        results.raw(i).kR_var  = kR_var_trajectory(:);


        results.raw(i).k_hist = processor.ekf_k.k_hist;
        results.raw(i).k_var  = processor.ekf_k.P_hist;

        tol = 1e-9;
        results.raw(i).at_bound_van    = abs(processor.ekf_van.current_params(:) - cfg.lb(:)) < tol ...
                                       | abs(processor.ekf_van.current_params(:) - cfg.ub(:)) < tol;
        results.raw(i).at_bound_gre    = abs(processor.ekf_gre.current_params(:) - cfg.lb(:)) < tol ...
                                       | abs(processor.ekf_gre.current_params(:) - cfg.ub(:)) < tol;
        results.raw(i).at_bound_loglin = abs(processor.ekf_loglin.current_params(:) - processor.ekf_loglin.lb(:)) < tol ...
                                       | abs(processor.ekf_loglin.current_params(:) - processor.ekf_loglin.ub(:)) < tol;
        results.raw(i).at_bound_2d     = abs(processor.ekf_2d_fim.current_params(:) - processor.ekf_2d_fim.lb(:)) < tol ...
                                       | abs(processor.ekf_2d_fim.current_params(:) - processor.ekf_2d_fim.ub(:)) < tol;
        results.raw(i).at_bound_1d     = abs(processor.ekf_k.k - processor.ekf_k.lb_k) < tol ...
                                       | abs(processor.ekf_k.k - processor.ekf_k.ub_k) < tol;
        
        results.metrics.van.MAE(i,1) = mae_v;
        results.metrics.gre.MAE(i,1) = mae_g;
        results.metrics.pop.MAE(i,1) = mae_p;
        results.metrics.kscale.MAE(i,1) = mae_k;
        results.metrics.loglin.MAE(i,1) = mae_ll;
        results.metrics.m2d.MAE(i,1) = mae_2d_fim;
        results.metrics.m2d_fim.MAE(i,1) = mae_2d_fim;
        
        results.raw(i).pers_van = processor.personalization.van;
        results.raw(i).pers_gre = processor.personalization.gre;
    catch ME
        fprintf('FAIL: %s\n', ME.message);
        results.fail_ids(end+1) = pid;
    end
end

if ~isempty(results.fail_ids)
    error('%d of %d cohort cases failed: %s', ...
        numel(results.fail_ids), numel(cohort), mat2str(results.fail_ids));
end
prov = provenance_stamp(cfg, DATA_FILE);
save(fullfile(results_dir, 'bis_analysis_results_v6_0.mat'), 'results', 'cfg', 'prov', '-v7.3');
%% =============== PARAMETER MULTIPLICITY ANALYSIS ==========
pid_test = 105;
eq_info = monte_carlo_equivalence_radius(results, pid_test, cfg);

fprintf('\n=== Equivalence analysis: patient %d ===\n', pid_test);
fprintf('Online theta = [%.2f, %.2f, %.2f, %.2f], MAE_ref = %.2f\n', ...
    eq_info.theta_ref, eq_info.mae_ref);

if ~eq_info.has_eq
    fprintf('No additional equivalent parameter sets found within +10%% MAE.\n');
else
    fprintf('Found %d random parameter sets with MAE <= 1.10 * MAE_ref.\n', eq_info.N_eq);
    fprintf('Equivalence radius (L2) mean = %.3f, max = %.3f\n', ...
        eq_info.radius_mean, eq_info.radius_max);
    fprintf('Example equivalent theta (first 5 rows):\n');
    disp(eq_info.theta_good(1:min(5,eq_info.N_eq), :));
end
%%
investigate_case(results, 278);
%%
plot_covariance_evolution(results, 101);
plot_fim_diagnostics(results, 101);
plot_kscale_comparison(results, 101);
plot_loglin_analysis(results, 101);
plot_2d_comparison(results, 101, cfg, results_dir);
plot_population_stats_boxplot(results);
%%
plot_identifiability_comparison(results, results_dir);
plot_1d_vs_2d(results, 105, results_dir);
plot_1d_vs_3d(results, 105, results_dir);
plot_parameter_clustering_boxplots(results, results_dir);
%%
compute_statistical_tests(results);
analyze_high_improvement_cases(results);
%% ===================== JOURNAL FIGURES =====================
fprintf('\n=== GENERATING JOURNAL FIGURES ===\n');

fig_dir = results_dir;

%% Figure 2: MAE Boxplots (all models)
generate_figure2_mae_boxplots(results, fig_dir);

%% Figure 3: MAE vs Model Dimension
generate_figure3_dimensionality(results, fig_dir);

%% Figure 4: Representative Case Study
generate_figure4_case_study(results, 107, cfg, fig_dir);

%% Figure 5: Monte Carlo Equivalence
generate_figure5_monte_carlo(eq_info, cfg, fig_dir);

% Figure 6: 2D Parameter Drift
generate_figure6_2d_drift(results, fig_dir);

%% Figure 7: FIM Diagnostics
%generate_figure7_fim_diagnostics(results, fig_dir);
generate_combined_identifiability_figure(results, cfg, fig_dir, 101);
%% Figure 8+: Case Studies and 4D Parameter Evolution
fprintf('Generating case studies and 4D parameter evolution figures...\n');
if exist('build_vec_from_results','file') || true
    vec = build_vec_from_results(results);
    cases = pick_case_studies(vec, results);
    plot_case_studies(cases, results, cfg, fig_dir);
    fprintf('  Saved: figure_cases_comparison.png\n');
    plot_case_study_params(cases, results, cfg, fig_dir);
    fprintf('  Saved: figure_4d_params.png\n');
end

%% ===================== JOURNAL TABLES =====================
fprintf('\n=== GENERATING JOURNAL TABLES ===\n');
generate_table1_demographics(results, fig_dir);
generate_table2_tuning(cfg, fig_dir);
generate_table3_accuracy(results, fig_dir);
generate_table4_equivalence(results, fig_dir);
generate_table5_equivalent_params(eq_info, pid_test, fig_dir);


fprintf('\n=== PERSONALIZATION SPEED ANALYSIS ===\n');
report_personalization_stats(results);
generate_simulation_summary(results, cfg);
prov = provenance_stamp(cfg, DATA_FILE);
save(fullfile(results_dir, 'bis_analysis_results_v6_0.mat'), 'results', 'cfg', 'eq_info', 'pid_test', 'prov', '-v7.3');
fprintf('\nV6.0 (Rigorous Estimator) Analysis complete. Results saved.\n');
