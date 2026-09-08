function res = run_outer(A, b, rho, mu, lambda, x0, K, Tref, Tsched_fun)
% RUN_OUTER Chay thuat toan diem gan ke tang toc khong chinh xac.

n = length(x0);
x = x0; u = x0;
Aparam = 1;     
eta1 = 0; delta1 = 0; delta2 = 0;

xs = zeros(n, K+1); xs(:,1) = x0;
Avec = zeros(K+1, 1); Avec(1) = Aparam;
eps1v = zeros(K, 1); eps2v = zeros(K, 1);
delta1v = zeros(K+1, 1); delta2v = zeros(K+1, 1);
Fvals = zeros(K+1, 1); Fvals(1) = F_obj(A, b, rho, mu, x0);
alphav = zeros(K,1);

conv_factor = sqrt((lambda*rho + 1)/(lambda*rho));

for k = 1:K
    % Tinh toan trong so alpha
    alpha = (sqrt((Aparam*lambda)^2 + 4*Aparam*lambda) - Aparam*lambda) / 2;
    y = (1 - alpha)*x + alpha*u;

    % Giai bai toan con
    Tk = min(max(Tsched_fun(k), 1), Tref);
    [xtrace, objtrace] = inner_fista_trace(A, b, rho, mu, lambda, y, x, Tref);

    xnew  = xtrace(:, Tk);
    dinner = max(objtrace(Tk) - min(objtrace), 0);

    % Cap nhat do chinh xac va sai so
    eps1 = sqrt(2*lambda*dinner);
    eps2 = eps1 * conv_factor;

    Anew = (1 - alpha) * Aparam;
    unew = u - (1/alpha) * (y - xnew);

    eta1   = eta1 + eps1/alpha;
    delta1 = (1 - alpha)*delta1 + (alpha*eta1)^2 / (2*lambda);
    delta2 = (1 - alpha)*delta2 + eps2^2 / (2*lambda);

    % Luu vet
    x = xnew; u = unew; Aparam = Anew;
    xs(:, k+1) = x; Avec(k+1) = Aparam;
    eps1v(k) = eps1; eps2v(k) = eps2;
    delta1v(k+1) = delta1; delta2v(k+1) = delta2;
    Fvals(k+1) = F_obj(A, b, rho, mu, x); alphav(k) = alpha;
end

res = struct('xs', xs, 'Avec', Avec, 'eps1', eps1v, 'eps2', eps2v, ...
    'delta1', delta1v, 'delta2', delta2v, 'Fvals', Fvals, 'alpha', alphav);
end
