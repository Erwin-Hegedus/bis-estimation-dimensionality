here = fileparts(mfilename('fullpath'));
addpath(genpath(here));
P = repro_paths();
cfg = P.cfg;
DATA_FILE = P.data_file;
pid = 105;

Q_scale = diag([1.0, 1.5, 0.4, 0.15].^2);

run_one = @(Qv) deal_one(cfg, Qv, DATA_FILE, pid);

[P0, x0, dxsat0] = run_one(0 * Q_scale);
[P1, x1, dxsat1] = run_one(1e3 * Q_scale);

fprintf('Q=0    : final P diag = [%.4g %.4g %.4g %.4g]\n', diag(P0));
fprintf('Q=1e3  : final P diag = [%.4g %.4g %.4g %.4g]\n', diag(P1));
fprintf('final params  Q=0   = [%.4f %.4f %.4f %.4f]\n', x0);
fprintf('final params  Q=1e3 = [%.4f %.4f %.4f %.4f]\n', x1);
fprintf('fraction of samples that updated:  Q=0 %.3f   Q=1e3 %.3f\n', dxsat0, dxsat1);
fprintf('P differs? %d   params differ? %d\n', any(abs(diag(P0)-diag(P1))>1e-9), any(abs(x0-x1)>1e-9));

function [Pfin, xfin, upd_frac] = deal_one(cfg, Qv, DATA_FILE, pid)
    [bis, remi_rate, prop_rate, time] = load_patient_record(pid, DATA_FILE);
    n = min([numel(bis), numel(remi_rate), numel(prop_rate), numel(time)]);
    bis = bis(1:n); prop_rate_raw = prop_rate(1:n); remi_rate_raw = remi_rate(1:n); time = time(1:n);
    demo = get_patient_demo_from_mat(pid, DATA_FILE);
    [pk_prop, pk_remi] = pk_from_demographics_schnider_minto(demo, cfg);
    [prop_rate, remi_rate, bis, time, ~, n, ~] = align_bis_to_infusion( ...
        bis, prop_rate_raw, remi_rate_raw, time, 0.3/60);
    cfg.Q_van = Qv;
    processor = init_processor(cfg, time(1), pk_prop, pk_remi, demo);
    nupd = 0;
    for k = 1:n
        prev_uc = processor.ekf_van.update_count;
        [~, ~, ~, ~, ~, processor] = process_sample(processor, time(k), bis(k), prop_rate(k), remi_rate(k));
        if processor.ekf_van.update_count > prev_uc
            nupd = nupd + 1;
        end
    end
    Pfin = processor.ekf_van.P;
    xfin = processor.ekf_van.current_params(:)';
    upd_frac = nupd / max(1,n);
end
