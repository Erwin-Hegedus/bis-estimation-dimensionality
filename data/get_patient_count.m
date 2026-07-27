function N = get_patient_count(matfile)
    S = load(matfile, 'patientDataFinal');
    if istable(S.patientDataFinal)
        N = height(S.patientDataFinal);
    else
        N = numel(S.patientDataFinal);
    end
end
