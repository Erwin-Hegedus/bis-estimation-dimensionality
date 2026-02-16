function N = get_patient_count_auto(matfile)
    S = load(matfile, 'patientDataFinal');
    if istable(S.patientDataFinal)
        N = height(S.patientDataFinal);
    else
        N = numel(S.patientDataFinal);
    end
end
