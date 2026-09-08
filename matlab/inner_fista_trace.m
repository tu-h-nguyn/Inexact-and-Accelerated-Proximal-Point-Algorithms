function [xtrace, objtrace] = inner_fista_trace(A, b, rho, mu, lambda, y, x0, Tmax)
% INNER_FISTA_TRACE Dung FISTA giai xap xi bai toan con.

L = norm(A)^2 + rho + 1/lambda;
n = length(x0);
z = x0; xk = x0; xkm1 = x0; t = 1;

xtrace   = zeros(n, Tmax);
objtrace = zeros(Tmax, 1);

for j = 1:Tmax
    % Buoc Gradient (Smooth)
    grad = A' * (A*z - b) + rho*z + (z - y)/lambda;

    % Buoc Proximal (non-smooth L1)
    xk = soft_threshold(z - grad/L, mu/L);

    % Buoc quan tinh (FISTA acceleration)
    tnext = (1 + sqrt(1 + 4*t^2)) / 2;
    z = xk + ((t - 1)/tnext) * (xk - xkm1);

    xkm1 = xk; t = tnext;

    % Luu vet
    xtrace(:, j) = xk;
    objtrace(j)  = 0.5*norm(A*xk - b)^2 + (rho/2)*norm(xk)^2 ...
        + mu*norm(xk, 1) + (1/(2*lambda))*norm(xk - y)^2;
end
end
