function h = updateVj(alpha, b, h0, Q)
% solving the QP problem.
%            min  1/2 v' * Q * v - b * v - sum(alpha .* log(v))
max_iter = 500;
rho = .995;
tol = .0001;

r = size(Q, 1);
s = alpha ./ h0;
Jk = zeros(2*r + 1, 2*r + 1);
Jk(1:r, 1:r) = Q;
Jk(1:r, r+1) = -1;
Jk(1:r,  r+2:end) = -eye(r);
Jk(r+1, 1:r) = 1;
Jk(r+2:end, 1:r) = diag(s);
Jk(r+2:end, r+2:end) = diag(h0);
vk = ones(2*r + 1, 1);
vk(1:r) = h0;
vk(r+2:end) = s;

for n_iter = 1:max_iter
    Fa = kkt_res(vk, alpha, b, Q);
    if max(abs(Fa)) < tol
        break;
    end
    d = linsolve(Jk, -Fa);
    step = vk ./ d;
    step(r+1) = -1;
    %step
    if all(step <= 0)
        % nofeasible direction
        break;
    else
        max_step = min(step(step > 0));
        step_length = min(1, max_step * rho);
    end
    vk = vk + step_length * d;
    assert(all(vk(1:r) > 0));
    Jk(r+2:end, r+2:end) = diag(vk(1:r));
    Jk(r+2:end, 1:r) = diag(vk(r+2:end));
end
h = vk(1:r);
end


function F_a = kkt_res(v, alpha, b, Q)
r = size(Q, 1);
F_a = zeros(2*r + 1, 1);
F_a(1:r) = Q * v(1:r) - v(r+1) - v(r+2:end) - b;
F_a(r+1) = sum(v(1:r)) - 1;
F_a(r+2:end) = v(1:r) .* v(r+2:end) - alpha;
end