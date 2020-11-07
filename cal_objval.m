function val = cal_objval(alpha, b, v, Q)
val = 0.5 * v' * Q * v - b' * v - sum(alpha .* log(v));
end