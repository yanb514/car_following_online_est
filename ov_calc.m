function [optimal_v] = ov_calc(a,b,hm,s)
    optimal_v = a*(tanh((s-hm)/b)+tanh(hm/b));
end

