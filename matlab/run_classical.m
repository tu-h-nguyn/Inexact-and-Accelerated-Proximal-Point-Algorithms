function Fvals = run_classical(A, b, rho, mu, lambda, x0, K, Tref)
% RUN_CLASSICAL Chay thuat toan PPA khong tang toc de so sanh.
n = length(x0); x = x0;
Fvals = zeros(K+1,1); Fvals(1) = F_obj(A, b, rho, mu, x0);

for k = 1:K
    [xtrace, ~] = inner_fista_trace(A, b, rho, mu, lambda, x, x, Tref);
    x = xtrace(:, Tref);
    Fvals(k+1) = F_obj(A, b, rho, mu, x);
end
end
