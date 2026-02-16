function pk = init_pk_state_online(p)
    pk = struct('params', p, 'x', zeros(3,1), 'Cp', 0, 'last_time', NaN);
    pk.k10 = p.CL / p.V1 / 60;
    pk.k12 = p.Q2 / p.V1 / 60;
    pk.k21 = p.Q2 / p.V2 / 60;
    pk.k13 = p.Q3 / p.V1 / 60;
    pk.k31 = p.Q3 / p.V3 / 60;
end
