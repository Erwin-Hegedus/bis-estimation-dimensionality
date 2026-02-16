function processor = update_induction_detection_simple(processor, bis_k, prop_rate_k, remi_rate_k)
    total_rate = prop_rate_k + remi_rate_k;
    switch processor.induction_status.phase
        case 'pre_induction'
            if total_rate > 1.0
                processor.induction_status.phase = 'induction';
                processor.induction_status.induction_start_sample = processor.sample_count;
            end
        case 'induction'
            if processor.sample_count - processor.induction_status.induction_start_sample > 300
                if bis_k < 60 && total_rate < 50
                    processor.induction_status.phase = 'maintenance';
                    processor.induction_status.induction_end_sample = processor.sample_count;
                end
            end
    end
end
