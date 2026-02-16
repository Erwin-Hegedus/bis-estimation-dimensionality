function demo = get_patient_demo_from_mat(index, matfile)
    S = load(matfile, 'patientDataFinal');
    T = S.patientDataFinal;
    demo = struct('Age', 50, 'Wt_kg', 70, 'Ht_cm', 170, 'Sex', 'M');
    
    if istable(T)
        if ismember('Age', T.Properties.VariableNames)
            demo.Age = double(T.Age(index));
        end
        if ismember('Weight_kg', T.Properties.VariableNames)
            demo.Wt_kg = double(T.Weight_kg(index));
        end
        if ismember('Height_cm', T.Properties.VariableNames)
            demo.Ht_cm = double(T.Height_cm(index));
        end
        if ismember('Sex', T.Properties.VariableNames)
            s = T.Sex(index);
            if iscell(s)
                s = s{1};
            end
            demo.Sex = string(s);
        end
    end
end
