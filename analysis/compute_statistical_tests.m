function stats_results = compute_statistical_tests(results)
% COMPUTE_STATISTICAL_TESTS - Comprehensive statistical analysis for journal
% Performs paired comparisons, confidence intervals, effect sizes, and corrections

    fprintf('\n');
    fprintf('================================================================================\n');
    fprintf('                    STATISTICAL ANALYSIS FOR JOURNAL                           \n');
    fprintf('================================================================================\n\n');
    
    % Extract MAE values
    mae_pop = results.metrics.pop.MAE(:);
    mae_1d = results.metrics.kscale.MAE(:);
    mae_2d = results.metrics.m2d.MAE(:);
    mae_2d_fim = results.metrics.m2d_fim.MAE(:);
    mae_3d = results.metrics.loglin.MAE(:);
    mae_4d = results.metrics.van.MAE(:);
    mae_4d_gre = results.metrics.gre.MAE(:);
    
    % Valid cases (all models have data)
    valid = ~isnan(mae_pop) & ~isnan(mae_1d) & ~isnan(mae_2d) & ...
            ~isnan(mae_2d_fim) & ~isnan(mae_3d) & ~isnan(mae_4d) & ~isnan(mae_4d_gre);
    
    N = sum(valid);
    fprintf('Valid cases for analysis: N = %d\n\n', N);
    
    if N < 10
        warning('Too few cases for reliable statistical inference (N=%d)', N);
    end
    
    % Filter to valid
    mae_pop = mae_pop(valid);
    mae_1d = mae_1d(valid);
    mae_2d = mae_2d(valid);
    mae_2d_fim = mae_2d_fim(valid);
    mae_3d = mae_3d(valid);
    mae_4d = mae_4d(valid);
    mae_4d_gre = mae_4d_gre(valid);
    
    % Store all MAEs in struct for easy access
    mae = struct('pop', mae_pop, 'm1d', mae_1d, 'm2d', mae_2d, ...
                 'm2d_fim', mae_2d_fim, 'm3d', mae_3d, 'm4d', mae_4d, 'm4d_gre', mae_4d_gre);
    
    %% ==================== DESCRIPTIVE STATISTICS ====================
    fprintf('--- DESCRIPTIVE STATISTICS (MAE in BIS units) ---\n\n');
    fprintf('%-15s %8s %8s %8s %12s %15s\n', 'Model', 'Mean', 'SD', 'Median', 'IQR', '95% CI');
    fprintf('%s\n', repmat('-', 1, 70));
    
    models = {'pop', 'm1d', 'm2d', 'm2d_fim', 'm3d', 'm4d', 'm4d_gre'};
    model_names = {'Population', '1D (k)', '2D (kP,kR)', '2D-FIM', '3D (LogLin)', '4D (Bouillon)', '4D (Greco)'};
    
    for i = 1:length(models)
        data = mae.(models{i});
        m = mean(data);
        s = std(data);
        med = median(data);
        iqr_val = [prctile(data, 25), prctile(data, 75)];
        ci = compute_ci(data, 0.95);
        
        fprintf('%-15s %8.2f %8.2f %8.2f %5.2f-%-5.2f %6.2f - %.2f\n', ...
            model_names{i}, m, s, med, iqr_val(1), iqr_val(2), ci(1), ci(2));
        
        stats_results.descriptive.(models{i}) = struct('mean', m, 'std', s, ...
            'median', med, 'iqr', iqr_val, 'ci95', ci, 'N', N);
    end
    fprintf('\n');
    
    %% ==================== PAIRWISE COMPARISONS ====================
    fprintf('--- PAIRWISE COMPARISONS (1D as reference) ---\n\n');
    
    % Define comparisons of interest
    comparisons = {
        'm1d', 'pop', '1D vs Pop';
        'm1d', 'm2d_fim', '1D vs 2D';
        'm1d', 'm3d', '1D vs 3D';
        'm1d', 'm4d', '1D vs 4D (Bouillon)';
        'm1d', 'm4d_gre', '1D vs 4D (Greco)';
        'm2d_fim', 'm3d', '2D vs 3D';
        'm2d_fim', 'm4d', '2D vs 4D';
        'm3d', 'm4d', '3D vs 4D';
        'm4d', 'm4d_gre', 'Bouillon vs Greco'
    };
    
    fprintf('%-20s %10s %10s %10s %10s %10s %10s\n', ...
        'Comparison', 'Mean Diff', '95% CI', 'p (Wilcox)', 'p (t-test)', 'Cohen d', 'Interpret');
    fprintf('%s\n', repmat('-', 1, 95));
    
    n_comparisons = size(comparisons, 1);
    p_values_wilcox = nan(n_comparisons, 1);
    p_values_ttest = nan(n_comparisons, 1);
    
    for i = 1:n_comparisons
        model_a = comparisons{i, 1};
        model_b = comparisons{i, 2};
        label = comparisons{i, 3};
        
        data_a = mae.(model_a);
        data_b = mae.(model_b);
        diff_ab = data_a - data_b;  % Positive = A worse than B
        
        % Mean difference and CI
        mean_diff = mean(diff_ab);
        ci_diff = compute_ci(diff_ab, 0.95);
        
        % Wilcoxon signed-rank test (non-parametric)
        [p_wilcox, ~] = signrank(data_a, data_b);
        p_values_wilcox(i) = p_wilcox;
        
        % Paired t-test (parametric)
        [~, p_ttest] = ttest(data_a, data_b);
        p_values_ttest(i) = p_ttest;
        
        % Effect size (Cohen's d for paired samples)
        cohens_d = mean_diff / std(diff_ab);
        effect_interp = interpret_cohens_d(cohens_d);
        
        fprintf('%-20s %+10.2f %+5.2f,%+5.2f %10.4f %10.4f %10.2f %10s\n', ...
            label, mean_diff, ci_diff(1), ci_diff(2), p_wilcox, p_ttest, cohens_d, effect_interp);
        
        stats_results.pairwise(i).comparison = label;
        stats_results.pairwise(i).mean_diff = mean_diff;
        stats_results.pairwise(i).ci_diff = ci_diff;
        stats_results.pairwise(i).p_wilcoxon = p_wilcox;
        stats_results.pairwise(i).p_ttest = p_ttest;
        stats_results.pairwise(i).cohens_d = cohens_d;
        stats_results.pairwise(i).effect_size = effect_interp;
    end
    
    fprintf('\n');
    
    %% ==================== MULTIPLE COMPARISON CORRECTION ====================
    fprintf('--- MULTIPLE COMPARISON CORRECTION ---\n\n');
    
    % Bonferroni correction
    alpha = 0.05;
    bonferroni_alpha = alpha / n_comparisons;
    bonferroni_sig = p_values_wilcox < bonferroni_alpha;
    
    % Benjamini-Hochberg FDR correction
    [p_fdr, fdr_sig] = fdr_correction(p_values_wilcox, alpha);
    
    fprintf('Original alpha: %.3f\n', alpha);
    fprintf('Bonferroni-corrected alpha: %.4f\n', bonferroni_alpha);
    fprintf('\n');
    
    fprintf('%-20s %12s %12s %12s %12s\n', 'Comparison', 'p (Wilcox)', 'Bonf. Sig', 'p (FDR)', 'FDR Sig');
    fprintf('%s\n', repmat('-', 1, 70));
    
    for i = 1:n_comparisons
        label = comparisons{i, 3};
        fprintf('%-20s %12.4f %12s %12.4f %12s\n', ...
            label, p_values_wilcox(i), bool2str(bonferroni_sig(i)), ...
            p_fdr(i), bool2str(fdr_sig(i)));
        
        stats_results.pairwise(i).bonferroni_sig = bonferroni_sig(i);
        stats_results.pairwise(i).fdr_p = p_fdr(i);
        stats_results.pairwise(i).fdr_sig = fdr_sig(i);
    end
    fprintf('\n');
    
    %% ==================== EQUIVALENCE TESTING (TOST) ====================
    fprintf('--- EQUIVALENCE TESTING (TOST) ---\n');
    fprintf('    H0: |MAE_1D - MAE_other| >= delta (clinically different)\n');
    fprintf('    H1: |MAE_1D - MAE_other| < delta (clinically equivalent)\n\n');
    
    delta = 1.0;  % 1 BIS unit equivalence margin
    fprintf('Equivalence margin (delta): %.1f BIS units\n\n', delta);
    
    fprintf('%-20s %12s %12s %12s\n', 'Comparison', 'Mean Diff', 'TOST p', 'Equivalent?');
    fprintf('%s\n', repmat('-', 1, 60));
    
    tost_comparisons = {'m2d_fim', 'm3d', 'm4d', 'm4d_gre'};
    tost_labels = {'1D vs 2D', '1D vs 3D', '1D vs 4D (Bou)', '1D vs 4D (Gre)'};
    
    for i = 1:length(tost_comparisons)
        data_other = mae.(tost_comparisons{i});
        diff_data = mae_1d - data_other;
        
        [p_tost, equivalent] = tost_paired(mae_1d, data_other, delta);
        mean_diff = mean(diff_data);
        
        fprintf('%-20s %+12.2f %12.4f %12s\n', ...
            tost_labels{i}, mean_diff, p_tost, bool2str(equivalent));
        
        stats_results.equivalence(i).comparison = tost_labels{i};
        stats_results.equivalence(i).mean_diff = mean_diff;
        stats_results.equivalence(i).tost_p = p_tost;
        stats_results.equivalence(i).equivalent = equivalent;
    end
    fprintf('\n');
    
    %% ==================== NON-INFERIORITY ANALYSIS ====================
    fprintf('--- NON-INFERIORITY ANALYSIS ---\n');
    fprintf('    1D is non-inferior if upper 95%% CI of (MAE_1D - MAE_4D) < margin\n\n');
    
    margin = 1.0;  % 1 BIS unit
    diff_1d_4d = mae_1d - mae_4d;
    ci_1d_4d = compute_ci(diff_1d_4d, 0.95);
    
    non_inferior = ci_1d_4d(2) < margin;
    
    fprintf('Non-inferiority margin: %.1f BIS units\n', margin);
    fprintf('MAE difference (1D - 4D): %.2f [95%% CI: %.2f, %.2f]\n', ...
        mean(diff_1d_4d), ci_1d_4d(1), ci_1d_4d(2));
    fprintf('Non-inferior: %s (upper CI %.2f < margin %.1f)\n', ...
        bool2str(non_inferior), ci_1d_4d(2), margin);
    fprintf('\n');
    
    stats_results.non_inferiority.margin = margin;
    stats_results.non_inferiority.mean_diff = mean(diff_1d_4d);
    stats_results.non_inferiority.ci = ci_1d_4d;
    stats_results.non_inferiority.is_non_inferior = non_inferior;
    
    %% ==================== WIN/LOSS ANALYSIS ====================
    fprintf('--- WIN/LOSS/TIE ANALYSIS (vs 4D Bouillon) ---\n\n');
    
    threshold = 0.5;  % Within 0.5 BIS = tie
    
    models_compare = {'m1d', 'm2d_fim', 'm3d'};
    models_compare_names = {'1D', '2D', '3D'};
    
    fprintf('%-10s %8s %8s %8s %12s\n', 'Model', 'Wins', 'Losses', 'Ties', 'Win Rate');
    fprintf('%s\n', repmat('-', 1, 50));
    
    for i = 1:length(models_compare)
        data_model = mae.(models_compare{i});
        diff_vs_4d = data_model - mae_4d;
        
        wins = sum(diff_vs_4d < -threshold);   % Model better by >threshold
        losses = sum(diff_vs_4d > threshold);  % Model worse by >threshold
        ties = N - wins - losses;
        win_rate = 100 * wins / N;
        
        fprintf('%-10s %8d %8d %8d %11.1f%%\n', ...
            models_compare_names{i}, wins, losses, ties, win_rate);
        
        stats_results.win_loss(i).model = models_compare_names{i};
        stats_results.win_loss(i).wins = wins;
        stats_results.win_loss(i).losses = losses;
        stats_results.win_loss(i).ties = ties;
        stats_results.win_loss(i).win_rate = win_rate;
    end
    fprintf('\n');
    
    %% ==================== FRIEDMAN TEST (Overall Comparison) ====================
    fprintf('--- OVERALL MODEL COMPARISON (Friedman Test) ---\n');
    fprintf('    Non-parametric alternative to repeated-measures ANOVA\n\n');
    
    % Combine all MAEs into matrix
    mae_matrix = [mae_pop, mae_1d, mae_2d_fim, mae_3d, mae_4d];
    
    try
        [p_friedman, tbl, stats_friedman] = friedman(mae_matrix, 1, 'off');
        
        fprintf('Friedman chi-square: %.2f\n', tbl{2,5});
        fprintf('Friedman p-value: %.4f\n', p_friedman);
        fprintf('Interpretation: %s\n', ...
            ternary(p_friedman < 0.05, 'Significant difference exists among models', ...
                                       'No significant overall difference'));
        
        stats_results.friedman.chisq = tbl{2,5};
        stats_results.friedman.p = p_friedman;
    catch ME
        fprintf('Friedman test failed: %s\n', ME.message);
        stats_results.friedman.chisq = NaN;
        stats_results.friedman.p = NaN;
    end
    fprintf('\n');
    
    %% ==================== SUMMARY FOR JOURNAL ====================
    fprintf('================================================================================\n');
    fprintf('                         JOURNAL-READY SUMMARY                                 \n');
    fprintf('================================================================================\n\n');
    
    fprintf('Sample size: N = %d patients\n\n', N);
    
    fprintf('Primary outcome (MAE, BIS units):\n');
    fprintf('  1D model: %.2f [95%% CI: %.2f-%.2f]\n', ...
        stats_results.descriptive.m1d.mean, ...
        stats_results.descriptive.m1d.ci95(1), ...
        stats_results.descriptive.m1d.ci95(2));
    fprintf('  4D model: %.2f [95%% CI: %.2f-%.2f]\n', ...
        stats_results.descriptive.m4d.mean, ...
        stats_results.descriptive.m4d.ci95(1), ...
        stats_results.descriptive.m4d.ci95(2));
    fprintf('\n');
    
    fprintf('Primary comparison (1D vs 4D):\n');
    idx_1d_4d = find(strcmp({stats_results.pairwise.comparison}, '1D vs 4D (Bouillon)'));
    fprintf('  Mean difference: %+.2f [95%% CI: %+.2f, %+.2f]\n', ...
        stats_results.pairwise(idx_1d_4d).mean_diff, ...
        stats_results.pairwise(idx_1d_4d).ci_diff(1), ...
        stats_results.pairwise(idx_1d_4d).ci_diff(2));
    fprintf('  Wilcoxon p-value: %.4f\n', stats_results.pairwise(idx_1d_4d).p_wilcoxon);
    fprintf('  Effect size (Cohen''s d): %.2f (%s)\n', ...
        stats_results.pairwise(idx_1d_4d).cohens_d, ...
        stats_results.pairwise(idx_1d_4d).effect_size);
    fprintf('\n');
    
    fprintf('Conclusion: ');
    if abs(stats_results.pairwise(idx_1d_4d).mean_diff) < 1.0 && ...
       stats_results.pairwise(idx_1d_4d).p_wilcoxon > 0.05
        fprintf('1D and 4D models show no clinically significant difference.\n');
    elseif stats_results.pairwise(idx_1d_4d).mean_diff > 0
        fprintf('4D model outperforms 1D model.\n');
    else
        fprintf('1D model outperforms 4D model.\n');
    end
    
    fprintf('\n================================================================================\n');
    
    stats_results.N = N;
end
