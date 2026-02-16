function art = init_artifact_gate_simple(~)
    art = struct('last_t', NaN, 'last_y', NaN, 'last_valid_y', NaN, 'bad_count', 0);
end
