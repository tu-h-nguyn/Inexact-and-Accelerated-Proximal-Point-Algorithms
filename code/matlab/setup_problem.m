function [A, b, x_true, par] = setup_problem(seed)
% SETUP_PROBLEM Sinh bai toan LASSO tong hop dung cho moi thi nghiem.
%
% Tach rieng khoi main_experiment.m de main_experiment.m va verify_bounds.m
% chac chan dung DUNG MOT bai toan, thay vi hai ban sao co the troi khoi nhau.
%
% LUU Y ve tai lap: MATLAB va Octave dung hai bo sinh so ngau nhien KHAC NHAU,
% nen cung mot `seed` van cho hai bai toan khac nhau giua hai moi truong. Vi
% the verify_bounds.m chi kiem tra cac khang dinh CAU TRUC (chan ly thuyet co
% bi vi pham khong, bac hoi tu co dung khong) -- nhung dieu phai dung tren moi
% thuc the cua bai toan -- chu khong so sanh tung chu so.

if nargin < 1 || isempty(seed)
    seed = 7;
end

% Octave khong co rng(); va randn('seed',...) chi gieo randn, khong gieo rand,
% trong khi randperm() lai dung rand. Phai gieo ca hai.
if exist('rng', 'builtin') == 5 || exist('rng', 'file') == 2
    try
        rng(seed);
    catch
        randn('seed', seed); rand('seed', seed);
    end
else
    randn('seed', seed); rand('seed', seed);
end

par = struct('n', 80, 'm', 30, 's', 8, ...
             'mu', 0.10, ...      % He so phat sparsity (L1)
             'rho', 0.02, ...     % He so loi manh (Ridge)
             'lambda', 2.0, ...   % Tham so prox ngoai
             'Tref', 200, ...     % So lap FISTA noi lam moc "gan dung"
             'K', 120, ...        % So buoc lap ngoai de ve do thi
             'Tref_ref', 500, ... % So lap noi cho moc tham chieu F*, x*
             'K_ref', 250, ...    % So buoc ngoai cho moc tham chieu
             'noise_floor', 1e-11);

A = randn(par.m, par.n) / sqrt(par.m);
x_true = zeros(par.n, 1);
idx = randperm(par.n, par.s);
x_true(idx) = 3 * randn(par.s, 1);
b = A * x_true + 0.01 * randn(par.m, 1);
end
