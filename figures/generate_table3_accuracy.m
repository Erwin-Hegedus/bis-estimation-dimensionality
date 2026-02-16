function generate_table3_accuracy(results, fig_dir)
% TABLE 3: Comparative Analysis of Prediction Accuracy

    fprintf('Generating Table 3: Accuracy Comparison...\n');
    
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
    
    fid = fopen(fullfile(fig_dir, 'table3_accuracy.txt'), 'w');
    fprintf(fid, 'TABLE 3: Comparative Analysis of Prediction Accuracy (N=%d)\n', N);
    fprintf(fid, '============================================================\n\n');
    fprintf(fid, 'Model              Dimension    Mean MAE (SD)      Parameters\n');
    fprintf(fid, '-----              ---------    -------------      ----------\n');
    fprintf(fid, 'Population         0D           %.2f (%.2f)        ---\n', mean(mae_pop), std(mae_pop));
    fprintf(fid, 'k-scale            1D           %.2f (%.2f)        k\n', mean(mae_k), std(mae_k));
    fprintf(fid, '(kP, kR)           2D           %.2f (%.2f)        kP, kR\n', mean(mae_2d), std(mae_2d));
    fprintf(fid, 'Log-linear         3D           %.2f (%.2f)        a0, aP, aR\n', mean(mae_ll), std(mae_ll));
    fprintf(fid, 'Bouillon           4D           %.2f (%.2f)        C50P, C50R, gamma, beta\n', mean(mae_van), std(mae_van));
    fprintf(fid, 'Greco              4D           %.2f (%.2f)        C50P, C50R, gamma, alpha\n', mean(mae_gre), std(mae_gre));
    fclose(fid);
    
    fprintf('  Saved: table3_accuracy.txt\n');
end
