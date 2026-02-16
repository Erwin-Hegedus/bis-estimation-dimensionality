function generate_table1_demographics(results, fig_dir)
% TABLE 1: Cohort Demographics and Clinical Characteristics

    fprintf('Generating Table 1: Demographics...\n');
    
    valid_idx = find(~isnan(results.metrics.van.MAE));
    N = length(valid_idx);
    
    ages = nan(N, 1);
    weights = nan(N, 1);
    heights = nan(N, 1);
    durations = nan(N, 1);
    sex_male = 0;
    
    for ii = 1:N
        idx = valid_idx(ii);
        if isfield(results, 'demographics') && idx <= length(results.demographics)
            ages(ii) = results.demographics(idx).Age;
            weights(ii) = results.demographics(idx).Weight;
            heights(ii) = results.demographics(idx).Height;
            durations(ii) = results.demographics(idx).Duration;
            if strcmpi(results.demographics(idx).Sex, 'M')
                sex_male = sex_male + 1;
            end
        end
    end
    
    % Compute statistics
    age_med = median(ages, 'omitnan');
    age_iqr = [prctile(ages, 25), prctile(ages, 75)];
    age_range = [min(ages), max(ages)];
    
    wt_med = median(weights, 'omitnan');
    wt_iqr = [prctile(weights, 25), prctile(weights, 75)];
    wt_range = [min(weights), max(weights)];
    
    ht_med = median(heights, 'omitnan');
    ht_iqr = [prctile(heights, 25), prctile(heights, 75)];
    ht_range = [min(heights), max(heights)];
    
    dur_med = median(durations, 'omitnan');
    dur_iqr = [prctile(durations, 25), prctile(durations, 75)];
    dur_range = [min(durations), max(durations)];
    
    pct_male = 100 * sex_male / N;
    
    % Write to file
    fid = fopen(fullfile(fig_dir, 'table1_demographics.txt'), 'w');
    fprintf(fid, 'TABLE 1: Cohort Demographics (N=%d valid cases)\n', N);
    fprintf(fid, '================================================\n\n');
    fprintf(fid, 'Characteristic          Median (IQR)           Range\n');
    fprintf(fid, '----------------        ----------------       --------\n');
    fprintf(fid, 'Age (years)             %.0f (%.0f-%.0f)            %.0f-%.0f\n', ...
            age_med, age_iqr(1), age_iqr(2), age_range(1), age_range(2));
    fprintf(fid, 'Weight (kg)             %.0f (%.0f-%.0f)            %.0f-%.0f\n', ...
            wt_med, wt_iqr(1), wt_iqr(2), wt_range(1), wt_range(2));
    fprintf(fid, 'Height (cm)             %.0f (%.0f-%.0f)          %.0f-%.0f\n', ...
            ht_med, ht_iqr(1), ht_iqr(2), ht_range(1), ht_range(2));
    fprintf(fid, 'Duration (min)          %.0f (%.0f-%.0f)           %.0f-%.0f\n', ...
            dur_med, dur_iqr(1), dur_iqr(2), dur_range(1), dur_range(2));
    fprintf(fid, 'Sex (M/F %%)             %.0f/%.0f\n', pct_male, 100-pct_male);
    fclose(fid);
    
    fprintf('  Saved: table1_demographics.txt\n');
end
