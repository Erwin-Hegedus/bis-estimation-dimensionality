function report_personalization_stats(results)
    fprintf('\n--- Personalization Detection Statistics ---\n');
    nP = numel(results.patient_id);
    valid_count = 0;
    van_times = nan(nP, 4);
    
    for ii = 1:nP
        if isfield(results.raw(ii), 'pers_van') && ~isempty(results.raw(ii).time)
            valid_count = valid_count + 1;
            d = results.raw(ii).pers_van.detected;
            t = results.raw(ii).pers_van.time;
            for j = 1:4
                if d(j)
                    van_times(ii, j) = t(j);
                end
            end
        end
    end
    
    fprintf('Analyzed %d valid cases (>30mins).\n', valid_count);
end
