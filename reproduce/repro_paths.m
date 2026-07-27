function P = repro_paths()
%REPRO_PATHS  Locate the inputs the reproduction scripts need.
%   P = REPRO_PATHS() returns a struct with the repository root, the patient
%   data file, the configuration used for the published run, the 209-case
%   cohort list and an output directory. Inputs that are not distributed with
%   the repository are reported by name here rather than failing later inside
%   a load.
%
%   Fields:
%     repo        repository root
%     data_file   300-case VitalDB extract, tracked at the repository root
%     cfg         configuration struct of the published run
%     cohort      209 patient indices of the analysed cohort
%     outdir      reproduce/output, created on first use

    here = fileparts(mfilename('fullpath'));
    P.repo = fileparts(here);

    P.outdir = fullfile(here, 'output');
    if ~exist(P.outdir, 'dir'), mkdir(P.outdir); end

    P.data_file = fullfile(P.repo, 'patientDataFinal_auto.mat');
    if ~exist(P.data_file, 'file')
        error('repro:missingData', ['patientDataFinal_auto.mat not found in %s.\n' ...
            'It is the 300-case VitalDB extract tracked at the repository root; ' ...
            'see reproduce/README.md.'], P.repo);
    end

    cfg_file = fullfile(here, 'cfg_canonical.mat');
    if ~exist(cfg_file, 'file')
        error('repro:missingCfg', 'cfg_canonical.mat not found in %s.', here);
    end
    S = load(cfg_file, 'cfg');
    P.cfg = S.cfg;

    csv = fullfile(here, 'cohort_caseids.csv');
    if ~exist(csv, 'file')
        error('repro:missingCohort', 'cohort_caseids.csv not found in %s.', here);
    end
    P.cohort = readmatrix(csv);
    P.cohort = P.cohort(:, 1);
end
