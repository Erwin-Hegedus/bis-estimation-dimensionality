function vec = build_vec_from_results_v3(results)
    vec.MAE_van = results.metrics.van.MAE(:);
    vec.MAE_gre = results.metrics.gre.MAE(:);
    vec.MAE_kscale = results.metrics.kscale.MAE(:);
    vec.MAE_loglin = results.metrics.loglin.MAE(:);
    vec.MAE_2d = results.metrics.m2d.MAE(:);
end
