function summary = generate_simulation_summary(results, cfg)
    fprintf('\n================================================================\n');
    fprintf('   IDENTIFIABILITY ANALYSIS: 0D vs 1D vs 2D vs 3D vs 4D        \n');
    fprintf('================================================================\n');
    
    valid_mask = ~isnan(results.metrics.van.MAE) & ...
                 ~isnan(results.metrics.pop.MAE) & ...
                 ~isnan(results.metrics.kscale.MAE) & ...
                 ~isnan(results.metrics.m2d.MAE) & ...
                 ~isnan(results.metrics.m2d_fim.MAE) & ...
                 ~isnan(results.metrics.loglin.MAE) & ...
                 ~isnan(results.metrics.gre.MAE);
    
    N_total = length(results.patient_id);
    N_valid = sum(valid_mask);
    
    mae_pop = results.metrics.pop.MAE(valid_mask);
    mae_van = results.metrics.van.MAE(valid_mask);
    mae_gre = results.metrics.gre.MAE(valid_mask);
    mae_k = results.metrics.kscale.MAE(valid_mask);
    mae_2d = results.metrics.m2d.MAE(valid_mask);
    mae_2d_fim = results.metrics.m2d_fim.MAE(valid_mask);
    mae_ll = results.metrics.loglin.MAE(valid_mask);
    
    fprintf('Processing %d valid cases (out of %d total).\n\n', N_valid, N_total);
    
    fprintf('=== MODEL PERFORMANCE COMPARISON ===\n\n');
    fprintf('| %-20s | %-6s | %-10s | %-10s |\n', 'Model', 'Dim', 'Mean MAE', 'Std MAE');
    fprintf('|----------------------|--------|------------|------------|\n');
    fprintf('| %-20s | %-6s | %10.2f | %10.2f |\n', 'Population (fixed)', '0D', mean(mae_pop), std(mae_pop));
    fprintf('| %-20s | %-6s | %10.2f | %10.2f |\n', 'k_scale (potency)', '1D', mean(mae_k), std(mae_k));
    fprintf('| %-20s | %-6s | %10.2f | %10.2f |\n', '2D (kP, kR)', '2D', mean(mae_2d), std(mae_2d));
    fprintf('| %-20s | %-6s | %10.2f | %10.2f |\n', '2D (kP, kR)', '2D', mean(mae_2d_fim), std(mae_2d_fim));
    fprintf('| %-20s | %-6s | %10.2f | %10.2f |\n', 'LogLin (a0,aP,aR)', '3D', mean(mae_ll), std(mae_ll));
    fprintf('| %-20s | %-6s | %10.2f | %10.2f |\n', 'Bouillon (full)', '4D', mean(mae_van), std(mae_van));
    fprintf('| %-20s | %-6s | %10.2f | %10.2f |\n', 'Greco (full)', '4D', mean(mae_gre), std(mae_gre));
    fprintf('\n');
    
    fprintf('=== PAIRWISE DIFFERENCES (mean ± std) ===\n\n');
    diff_pop_k = mae_pop - mae_k;
    diff_pop_2d = mae_pop - mae_2d;
    diff_k_2d = mae_k - mae_2d;
    diff_k_van = mae_k - mae_van;
    diff_2d_van = mae_2d - mae_van;
    diff_ll_van = mae_ll - mae_van;
    diff_van_gre = mae_van - mae_gre;
    
    fprintf('Pop - K1D:    %+.2f ± %.2f (K1D improvement over Pop)\n', mean(diff_pop_k), std(diff_pop_k));
    fprintf('Pop - 2D:     %+.2f ± %.2f (2D improvement over Pop)\n', mean(diff_pop_2d), std(diff_pop_2d));
    fprintf('K1D - 2D:     %+.2f ± %.2f (2D vs K1D)\n', mean(diff_k_2d), std(diff_k_2d));
    fprintf('K1D - Van4D:  %+.2f ± %.2f (Van4D vs K1D)\n', mean(diff_k_van), std(diff_k_van));
    fprintf('2D - Van4D:   %+.2f ± %.2f (Van4D vs 2D)\n', mean(diff_2d_van), std(diff_2d_van));
    fprintf('LL3D - Van4D: %+.2f ± %.2f (Van4D vs LL3D)\n', mean(diff_ll_van), std(diff_ll_van));
    fprintf('Van4D - Gre4D: %+.2f ± %.2f (Greco vs Bouillon)\n', mean(diff_van_gre), std(diff_van_gre));
    fprintf('\n');
    
    fprintf('=== IDENTIFIABILITY CONCLUSION ===\n\n');
    
    thresh = 1.0;
    
    n_k_eq_van = sum(abs(diff_k_van) <= thresh);
    n_2d_eq_van = sum(abs(diff_2d_van) <= thresh);
    n_k_eq_2d = sum(abs(diff_k_2d) <= thresh);
    
    fprintf('Cases where K(1D) ≈ Van(4D) (within %.1f BIS): %d/%d (%.0f%%)\n', ...
        thresh, n_k_eq_van, N_valid, 100*n_k_eq_van/N_valid);
    fprintf('Cases where 2D ≈ Van(4D) (within %.1f BIS): %d/%d (%.0f%%)\n', ...
        thresh, n_2d_eq_van, N_valid, 100*n_2d_eq_van/N_valid);
    fprintf('Cases where K(1D) ≈ 2D (within %.1f BIS): %d/%d (%.0f%%)\n', ...
        thresh, n_k_eq_2d, N_valid, 100*n_k_eq_2d/N_valid);
    fprintf('\n');
    
    if mean(abs(diff_k_2d)) < 0.5 && mean(abs(diff_k_van)) < 1.5
        fprintf('╔════════════════════════════════════════════════════════════════╗\n');
        fprintf('║  CONCLUSION: PD parameters are NOT individually identifiable!  ║\n');
        fprintf('║                                                                ║\n');
        fprintf('║  Evidence:                                                     ║\n');
        fprintf('║  - 1D, 2D, 3D, and 4D models achieve similar MAE               ║\n');
        fprintf('║  - Splitting k into kP and kR provides NO significant gain     ║\n');
        fprintf('║  - The BIS+TCI structure supports only ~1 effective DoF        ║\n');
        fprintf('║  - Additional parameters do not improve prediction accuracy   ║\n');
        fprintf('╚════════════════════════════════════════════════════════════════╝\n');
    elseif mean(diff_k_2d) > 0.5
        fprintf('╔════════════════════════════════════════════════════════════════╗\n');
        fprintf('║  FINDING: 2D model (kP, kR) provides improvement over 1D       ║\n');
        fprintf('║                                                                ║\n');
        fprintf('║  Mean improvement: %.2f BIS units                              ║\n', mean(diff_k_2d));
        fprintf('║  This suggests propofol/remi sensitivity may be independently  ║\n');
        fprintf('║  identifiable under certain excitation conditions.             ║\n');
        fprintf('╚════════════════════════════════════════════════════════════════╝\n');
    end
    
    summary.mae_pop = mean(mae_pop);
    summary.mae_k = mean(mae_k);
    summary.mae_2d = mean(mae_2d);
    summary.mae_ll = mean(mae_ll);
    summary.mae_van = mean(mae_van);
    summary.mae_gre = mean(mae_gre);
    summary.N_valid = N_valid;
end
