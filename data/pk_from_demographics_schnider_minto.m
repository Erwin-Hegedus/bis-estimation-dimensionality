function [pk_prop, pk_remi] = pk_from_demographics_schnider_minto(demo, cfg)
    age = demo.Age;
    wt = demo.Wt_kg;
    ht = demo.Ht_cm;
    sex_male = strncmpi(demo.Sex, 'M', 1);
    
    if sex_male
        LBM = 1.1 * wt - 128 * (wt / ht)^2;
    else
        LBM = 1.07 * wt - 148 * (wt / ht)^2;
    end
    LBM = max(30, min(wt, LBM));
    
    pk_prop = cfg.pk.prop;
    V1 = 4.27;
    V2 = 18.9 - 0.391 * (age - 53);
    V3 = 238;
    CL = 1.89 + 0.0456 * (wt - 77) - 0.0681 * (LBM - 59) + 0.0264 * (ht - 177);
    Q2 = 1.29 - 0.024 * (age - 53);
    Q3 = 0.836;
    pk_prop.V1 = V1;
    pk_prop.V2 = max(V2, 5);
    pk_prop.V3 = V3;
    pk_prop.CL = max(CL, 0.5);
    pk_prop.Q2 = max(Q2, 0.2);
    pk_prop.Q3 = Q3;
    
    pk_remi = cfg.pk.remi;
    V1r = 5.1;
    V2r = 9.82 - 0.0811 * (age - 40);
    V3r = 5.42;
    CLr = 2.6 - 0.0162 * (age - 40) + 0.0191 * (LBM - 55);
    Q2r = 2.05 - 0.0301 * (age - 40);
    Q3r = 0.076;
    pk_remi.V1 = V1r;
    pk_remi.V2 = max(V2r, 3);
    pk_remi.V3 = V3r;
    pk_remi.CL = max(CLr, 0.5);
    pk_remi.Q2 = max(Q2r, 0.3);
    pk_remi.Q3 = Q3r;
end
