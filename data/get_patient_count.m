function N = get_patient_count(matfile)
    persistent T cached_file;
    if isempty(T) || ~strcmp(cached_file, matfile)
        S = load(matfile, 'patientDataFinal');
        T = S.patientDataFinal;
        cached_file = matfile;
    end
    if istable(T)
        N = height(T);
    else
        N = numel(T);
    end
end
