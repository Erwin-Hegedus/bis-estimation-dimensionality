function improvement_stats = analyze_high_improvement_cases(results)
% ANALYZE_HIGH_IMPROVEMENT_CASES - Find cases where personalization helps most
% Identifies patients with high population MAE and significant improvement

    fprintf('\n');
    fprintf('================================================================================\n');
    fprintf('              HIGH IMPROVEMENT CASES ANALYSIS                                  \n');
    fprintf('================================================================================\n\n');

    mae_pop = results.metrics.pop.MAE(:);
    mae_1d = results.metrics.kscale.MAE(:);
    mae_2d = results.metrics.m2d_fim.MAE(:);
    mae_4d = results.metrics.van.MAE(:);
    
    valid = ~isnan(mae_pop) & ~isnan(mae_1d) & ~isnan(mae_2d) & ~isnan(mae_4d);
    
    idx_valid = find(valid);
    N = sum(valid);
    
    mae_pop = mae_pop(valid);
    mae_1d = mae_1d(valid);
    mae_2d = mae_2d(valid);
    mae_4d = mae_4d(valid);
    
    % Improvement (positive = model is better than population)
    improv_1d = mae_pop - mae_1d;
    improv_2d = mae_pop - mae_2d;
    improv_4d = mae_pop - mae_4d;
    
    fprintf('--- OVERALL IMPROVEMENT STATISTICS ---\n\n');
    fprintf('%-20s %10s %10s %10s\n', 'Metric', '1D', '2D', '4D');
    fprintf('%s\n', repmat('-', 1, 55));
    fprintf('%-20s %10.2f %10.2f %10.2f\n', 'Mean improvement', mean(improv_1d), mean(improv_2d), mean(improv_4d));
    fprintf('%-20s %10.2f %10.2f %10.2f\n', 'Max improvement', max(improv_1d), max(improv_2d), max(improv_4d));
    fprintf('%-20s %10d %10d %10d\n', 'Cases improved', sum(improv_1d > 0), sum(improv_2d > 0), sum(improv_4d > 0));
    fprintf('%-20s %10.1f%% %9.1f%% %9.1f%%\n', 'Percent improved', 100*mean(improv_1d > 0), 100*mean(improv_2d > 0), 100*mean(improv_4d > 0));
    fprintf('\n');
    
    %% High population MAE cases
    thresholds = [8, 10, 12, 15];
    
    fprintf('--- CASES WITH HIGH POPULATION MAE ---\n\n');
    fprintf('%-12s %8s %12s %12s %12s\n', 'Pop MAE >=', 'N cases', '1D improv', '2D improv', '4D improv');
    fprintf('%s\n', repmat('-', 1, 60));
    
    for thresh = thresholds
        mask = mae_pop >= thresh;
        n_cases = sum(mask);
        
        if n_cases > 0
            mean_imp_1d = mean(improv_1d(mask));
            mean_imp_2d = mean(improv_2d(mask));
            mean_imp_4d = mean(improv_4d(mask));
            fprintf('%-12d %8d %12.2f %12.2f %12.2f\n', thresh, n_cases, mean_imp_1d, mean_imp_2d, mean_imp_4d);
        else
            fprintf('%-12d %8d %12s %12s %12s\n', thresh, 0, '-', '-', '-');
        end
    end
    fprintf('\n');
    
    %% Detailed analysis for Pop MAE >= 10
    thresh_main = 10;
    mask_high = mae_pop >= thresh_main;
    n_high = sum(mask_high);
    
    fprintf('--- DETAILED: CASES WITH POPULATION MAE >= %.0f (N=%d) ---\n\n', thresh_main, n_high);
    
    if n_high > 0
        fprintf('Population MAE:  %.2f +/- %.2f [range: %.1f - %.1f]\n', ...
            mean(mae_pop(mask_high)), std(mae_pop(mask_high)), ...
            min(mae_pop(mask_high)), max(mae_pop(mask_high)));
        fprintf('1D MAE:          %.2f +/- %.2f [range: %.1f - %.1f]\n', ...
            mean(mae_1d(mask_high)), std(mae_1d(mask_high)), ...
            min(mae_1d(mask_high)), max(mae_1d(mask_high)));
        fprintf('2D MAE:          %.2f +/- %.2f [range: %.1f - %.1f]\n', ...
            mean(mae_2d(mask_high)), std(mae_2d(mask_high)), ...
            min(mae_2d(mask_high)), max(mae_2d(mask_high)));
        fprintf('4D MAE:          %.2f +/- %.2f [range: %.1f - %.1f]\n', ...
            mean(mae_4d(mask_high)), std(mae_4d(mask_high)), ...
            min(mae_4d(mask_high)), max(mae_4d(mask_high)));
        fprintf('\n');
        
        fprintf('Improvement over population:\n');
        fprintf('  1D: %.2f +/- %.2f BIS (%.0f%% reduction)\n', ...
            mean(improv_1d(mask_high)), std(improv_1d(mask_high)), ...
            100 * mean(improv_1d(mask_high)) / mean(mae_pop(mask_high)));
        fprintf('  2D: %.2f +/- %.2f BIS (%.0f%% reduction)\n', ...
            mean(improv_2d(mask_high)), std(improv_2d(mask_high)), ...
            100 * mean(improv_2d(mask_high)) / mean(mae_pop(mask_high)));
        fprintf('  4D: %.2f +/- %.2f BIS (%.0f%% reduction)\n', ...
            mean(improv_4d(mask_high)), std(improv_4d(mask_high)), ...
            100 * mean(improv_4d(mask_high)) / mean(mae_pop(mask_high)));
        fprintf('\n');
        
        % Cases where 1D achieves clinically meaningful improvement
        meaningful_thresh = 3;  % 3 BIS improvement
        n_meaningful_1d = sum(improv_1d(mask_high) >= meaningful_thresh);
        n_meaningful_2d = sum(improv_2d(mask_high) >= meaningful_thresh);
        n_meaningful_4d = sum(improv_4d(mask_high) >= meaningful_thresh);
        
        fprintf('Cases with >= %.0f BIS improvement:\n', meaningful_thresh);
        fprintf('  1D: %d/%d (%.0f%%)\n', n_meaningful_1d, n_high, 100*n_meaningful_1d/n_high);
        fprintf('  2D: %d/%d (%.0f%%)\n', n_meaningful_2d, n_high, 100*n_meaningful_2d/n_high);
        fprintf('  4D: %d/%d (%.0f%%)\n', n_meaningful_4d, n_high, 100*n_meaningful_4d/n_high);
        fprintf('\n');
    end
    
    %% Relative improvement analysis
    fprintf('--- RELATIVE IMPROVEMENT (stratified by population MAE) ---\n\n');
    
    bins = [0, 6, 8, 10, 15, 100];
    bin_labels = {'<6', '6-8', '8-10', '10-15', '>15'};
    
    fprintf('%-12s %8s %15s %15s %15s\n', 'Pop MAE', 'N', '1D improvement', '2D improvement', '4D improvement');
    fprintf('%s\n', repmat('-', 1, 70));
    
    for i = 1:length(bins)-1
        mask = mae_pop >= bins(i) & mae_pop < bins(i+1);
        n_bin = sum(mask);
        
        if n_bin > 0
            fprintf('%-12s %8d %10.2f +/- %-4.1f %10.2f +/- %-4.1f %10.2f +/- %-4.1f\n', ...
                bin_labels{i}, n_bin, ...
                mean(improv_1d(mask)), std(improv_1d(mask)), ...
                mean(improv_2d(mask)), std(improv_2d(mask)), ...
                mean(improv_4d(mask)), std(improv_4d(mask)));
        end
    end
    fprintf('\n');
    
    %% Best individual cases
    fprintf('--- TOP 10 CASES WITH HIGHEST 1D IMPROVEMENT ---\n\n');
    [~, sort_idx] = sort(improv_1d, 'descend');
    top_n = min(10, N);
    
    fprintf('%-10s %10s %10s %10s %10s %10s\n', 'Patient', 'Pop MAE', '1D MAE', '2D MAE', '4D MAE', 'Improv');
    fprintf('%s\n', repmat('-', 1, 65));
    
    for i = 1:top_n
        j = sort_idx(i);
        pid = results.patient_id(idx_valid(j));
        fprintf('%-10d %10.2f %10.2f %10.2f %10.2f %10.2f\n', ...
            pid, mae_pop(j), mae_1d(j), mae_2d(j), mae_4d(j), improv_1d(j));
    end
    fprintf('\n');
    
    %% Summary for paper
    fprintf('================================================================================\n');
    fprintf('                         JOURNAL-READY SUMMARY                                 \n');
    fprintf('================================================================================\n\n');
    
    n_pop_high = sum(mae_pop >= 10);
    pct_pop_high = 100 * n_pop_high / N;
    
    if n_pop_high > 0
        mean_improv_1d_high = mean(improv_1d(mae_pop >= 10));
        pct_reduction = 100 * mean_improv_1d_high / mean(mae_pop(mae_pop >= 10));
        
        fprintf('Of %d patients, %d (%.0f%%) had population MAE >= 10 BIS.\n', N, n_pop_high, pct_pop_high);
        fprintf('In these poorly-fit patients, the 1D model reduced MAE by %.1f BIS (%.0f%% reduction).\n', ...
            mean_improv_1d_high, pct_reduction);
        fprintf('\n');
    end
    
    % Cases where personalization is essential
    n_essential = sum(mae_pop >= 10 & mae_1d < 7);
    fprintf('Patients requiring personalization (Pop MAE >= 10, reduced to < 7 by 1D): %d (%.0f%%)\n', ...
        n_essential, 100*n_essential/N);
    
    fprintf('================================================================================\n\n');
    
    % Store results
    improvement_stats.N = N;
    improvement_stats.mae_pop = mae_pop;
    improvement_stats.mae_1d = mae_1d;
    improvement_stats.mae_2d = mae_2d;
    improvement_stats.mae_4d = mae_4d;
    improvement_stats.improv_1d = improv_1d;
    improvement_stats.improv_2d = improv_2d;
    improvement_stats.improv_4d = improv_4d;
    improvement_stats.n_pop_high = n_pop_high;
    improvement_stats.idx_valid = idx_valid;
end
