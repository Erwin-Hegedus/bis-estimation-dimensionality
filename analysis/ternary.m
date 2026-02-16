function result = ternary(condition, true_val, false_val)
% Ternary operator
    if condition
        result = true_val;
    else
        result = false_val;
    end
end