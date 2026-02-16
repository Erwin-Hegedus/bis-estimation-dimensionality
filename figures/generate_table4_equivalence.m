function generate_table4_equivalence(results, fig_dir)
% TABLE 4: Model Equivalence Rates and Pairwise Differences

    fprintf('Generating Table 4: Equivalence Rates...\n');
    
    valid_mask = ~isnan(results.metrics.van.MAE) & ...
                 ~isnan(results.metrics.gre.MAE) & ...
                 ~isnan(results.metrics.pop.MAE) & ...
                 ~isnan(results.metrics.kscale.MAE) & ...
                 ~isnan(results.metrics.m2d.MAE) & ...
                 ~isnan(results.metrics.loglin.MAE);
    
    N = sum(valid_mask);
    
    mae_pop = results.metrics.pop.MAE(valid_mask);
    mae_k = results.metrics.kscale.MAE(valid_mask);
    mae_2d = results.metrics.m2d.MAE(valid_mask);
    mae_ll = results.metrics.loglin.MAE(valid_mask);
    mae_van = results.metrics.van.MAE(valid_mask);
    mae_gre = results.metrics.gre.MAE(valid_mask);
    
    % Compute pairwise differences
    diff_1d_4d = mae_k - mae_van;
    diff_2d_4d = mae_2d - mae_van;
    diff_1d_2d = mae_k - mae_2d;
    diff_3d_4d = mae_ll - mae_van;
    diff_van_gre = mae_van - mae_gre;
    
    % Compute equivalence rates (within 1 BIS)
    thresh = 1.0;
    eq_1d_4d = sum(abs(diff_1d_4d) <= thresh);
    eq_2d_4d = sum(abs(diff_2d_4d) <= thresh);
    eq_1d_2d = sum(abs(diff_1d_2d) <= thresh);
    eq_3d_4d = sum(abs(diff_3d_4d) <= thresh);
    eq_van_gre = sum(abs(diff_van_gre) <= thresh);
    
    fid = fopen(fullfile(fig_dir, 'table4_equivalence.txt'), 'w');
    fprintf(fid, 'TABLE 4: Model Equivalence Rates and Pairwise Differences (N=%d)\n', N);
    fprintf(fid, '================================================================\n\n');
    fprintf(fid, 'A) Equivalence Rates (|Delta MAE| <= 1 BIS)\n');
    fprintf(fid, '-------------------------------------------\n');
    fprintf(fid, '1D vs 4D Bouillon:    %d/%d (%.0f%%)\n', eq_1d_4d, N, 100*eq_1d_4d/N);
    fprintf(fid, '2D vs 4D Bouillon:    %d/%d (%.0f%%)\n', eq_2d_4d, N, 100*eq_2d_4d/N);
    fprintf(fid, '1D vs 2D:             %d/%d (%.0f%%)\n', eq_1d_2d, N, 100*eq_1d_2d/N);
    fprintf(fid, '3D vs 4D Bouillon:    %d/%d (%.0f%%)\n', eq_3d_4d, N, 100*eq_3d_4d/N);
    fprintf(fid, 'Bouillon vs Greco:    %d/%d (%.0f%%)\n', eq_van_gre, N, 100*eq_van_gre/N);
    fprintf(fid, '\n');
    fprintf(fid, 'B) Pairwise MAE Differences (mean +/- SD)\n');
    fprintf(fid, '-----------------------------------------\n');
    fprintf(fid, '1D - 4D Bouillon:     %+.2f +/- %.2f\n', mean(diff_1d_4d), std(diff_1d_4d));
    fprintf(fid, '2D - 4D Bouillon:     %+.2f +/- %.2f\n', mean(diff_2d_4d), std(diff_2d_4d));
    fprintf(fid, '1D - 2D:              %+.2f +/- %.2f\n', mean(diff_1d_2d), std(diff_1d_2d));
    fprintf(fid, '3D - 4D Bouillon:     %+.2f +/- %.2f\n', mean(diff_3d_4d), std(diff_3d_4d));
    fprintf(fid, 'Bouillon - Greco:     %+.2f +/- %.2f\n', mean(diff_van_gre), std(diff_van_gre));
    fclose(fid);
    
    fprintf('  Saved: table4_equivalence.txt\n');
end
