% MAIN_EXPERIMENT
% Minh hoa so cho thuat toan diem gan ke khong chinh xac va tang toc (IAPPA).
%
% Bai toan: elastic-net / LASSO tong hop  F(x) = 1/2||Ax-b||^2 + rho/2||x||^2 + mu||x||_1
% Bai toan con prox duoc giai xap xi bang FISTA noi; so vong lap noi T_k dong vai tro
% "nut van" sinh sai so co kiem soat. Script xuat 3 hinh vao thu muc figures/ o goc repo.
%
% Cach chay:  >> cd code/matlab; main_experiment          (MATLAB)
%             $ cd code/matlab && octave --no-gui --quiet --eval "main_experiment"
%
% Chay duoc tren ca MATLAB lan GNU Octave: export_fig_pdf.m tu chon
% exportgraphics hay print tuy moi truong, nen khong can giay phep MATLAB.

clear; close all; clc;

% Ve hinh duoi Octave khong phai luc nao cung chay duoc: ft_text_renderer cua
% no hong tren mot so he thong ("unable to load font .../FreeSans.otf" roi
% "invalid bounding box"). Da thu ca gnuplot-nox lan chi dinh font DejaVu tren
% hai moi truong khac nhau, deu khong qua. PHAN TINH TOAN thi chay binh thuong.
%
% Nen: tinh toan luon chay va luon in ket qua; phan ve hinh duoc bao trong
% try/catch de Octave khong che mat ket qua so bang mot loi do hoa.
CAN_PLOT = true;
if exist('OCTAVE_VERSION', 'builtin')
    try
        h = figure('visible', 'off'); plot(1:2); close(h);
    catch
        CAN_PLOT = false;
        warning(['Octave khong ve duoc hinh tren may nay; bo qua phan do thi. ' ...
                 'Ket qua so o tren van day du. Dung MATLAB de sinh lai ' ...
                 'figures/*.pdf cho bao cao.']);
    end
end

% Thu muc hinh: <goc repo>/figures  (script nam tai <goc repo>/code/matlab)
here   = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
figdir = fullfile(here, '..', '..', 'figures');
if ~exist(figdir, 'dir'), mkdir(figdir); end

%% 1) Khoi tao va sinh du lieu bai toan
% Bai toan lay tu setup_problem.m -- CHUNG voi verify_bounds.m, de hinh ve va
% cong kiem tra khong the noi ve hai bai toan khac nhau. setup_problem cung
% gieo CA rand lan randn: randperm() dung rand, ma randn('seed',...) khong
% gieo rand, nen ban truoc day khong tai lap duoc tren Octave.
[A, b, ~, par] = setup_problem(7);
n = par.n; m = par.m; s = par.s;
mu = par.mu; rho = par.rho; lambda = par.lambda;
Tref = par.Tref; K = par.K;
Tref_ref = par.Tref_ref; K_ref = par.K_ref;
NOISE_FLOOR = par.noise_floor;
x0 = zeros(n, 1);

fprintf('Kich thuoc bai toan: n=%d, m=%d, do thua s=%d\n', n, m, s);

%% 2) Uoc luong F*, x* (chay moc tham chieu do chinh xac cao)
sched_exact_ref = @(k) Tref_ref;
res_ref = run_outer(A, b, rho, mu, lambda, x0, K_ref, Tref_ref, sched_exact_ref);
F_star = res_ref.Fvals(end);
x_star = res_ref.xs(:, end);
fprintf('F* (tham chieu, %d buoc ngoai x %d lap noi): %.10e\n', K_ref, Tref_ref, F_star);

%% 3) Thu nghiem cac lich trinh lap noi (mo phong sai so)
sched = schedules(Tref);

names = fieldnames(sched);
R = struct();
for i = 1:numel(names)
    nm = names{i};
    R.(nm) = run_outer(A, b, rho, mu, lambda, x0, K, Tref, sched.(nm));
    fprintf('%s: F(x_K)-F* = %.4e\n', nm, R.(nm).Fvals(end) - F_star);
end

%% 4) Duong co dien (khong tang toc)
Fvals_classical = run_classical(A, b, rho, mu, lambda, x0, K, Tref);
fprintf('Co dien (khong tang toc): F(x_K)-F* = %.4e\n', Fvals_classical(end) - F_star);

%% 5) Kiem chung toan hoc cac chan ly thuyet
A0 = R.S4.Avec(1);
phi0_gap = (F_obj(A,b,rho,mu,x0) - F_star) + (A0/2)*norm(x_star - x0)^2;

k_axis   = (0:K)';
beta_k   = R.S4.Avec / R.S4.Avec(1);
bound1   = beta_k .* phi0_gap + R.S4.delta1;
bound2   = beta_k .* phi0_gap + R.S4.delta2;
actual_S4 = R.S4.Fvals - F_star;

viol1 = sum(actual_S4 > bound1 + 1e-9);
viol2 = sum(actual_S4 > bound2 + 1e-9);
fprintf('Vi pham chan Loai 1 (S4): %d / %d\n', viol1, K+1);
fprintf('Vi pham chan Loai 2 (S4): %d / %d\n', viol2, K+1);

if CAN_PLOT
    %% 6) Truc quan hoa va danh gia thuc nghiem
    clipv = @(v) max(v, NOISE_FLOOR);
    kk = (1:K)';

    % -- HINH A: So sanh su hoi tu --
    fig1 = figure('Position',[100 100 640 480]);
    loglog(kk, clipv(Fvals_classical(2:end)-F_star), '-', 'LineWidth',1.6, 'Color',[0.3 0.3 0.3]); hold on;
    loglog(kk, clipv(R.S1.Fvals(2:end)-F_star), '-', 'LineWidth',1.7, 'Color',[0.64 0.08 0.18]);
    loglog(kk, clipv(R.S2.Fvals(2:end)-F_star), '-.', 'LineWidth',1.6, 'Color',[0.85 0.33 0.10]);
    loglog(kk, clipv(R.S4.Fvals(2:end)-F_star), '--', 'LineWidth',1.6, 'Color',[0.47 0.67 0.19]);
    loglog(kk, clipv(R.S5.Fvals(2:end)-F_star), '-', 'LineWidth',1.8, 'Color',[0 0.45 0.74]);
    loglog(kk, 0.5*(F_obj(A,b,rho,mu,x0)-F_star)./kk, ':', 'LineWidth', 1.2, 'Color', [0 0 0]);
    loglog(kk, 2*(F_obj(A,b,rho,mu,x0)-F_star)./kk.^2, ':', 'LineWidth', 1.2, 'Color', [0.5 0.5 0.5]);
    hold off; grid on; ylim([NOISE_FLOOR/2, 1e2]);
    xlabel('Buoc ngoai k'); ylabel('F(x_k) - F^*');
    legend('PPA Co dien', 'T_k=1 (Cuc tho)', 'T_k=4 (Tho)', 'T_k ~ sqrt(k)', 'T_k = T_{ref} (Chinh xac)', ...
        'Tham chieu O(1/k)', 'Tham chieu O(1/k^2)', 'Location','southwest');
    title('So sanh toc do hoi tu');

    filename1 = fullfile(figdir, 'convergence-comparison.pdf');
    export_fig_pdf(fig1, filename1);

    %% ================= HINH B: chan ly thuyet vs sai so thuc =================
    fig2 = figure('Position',[100 100 640 480]);
    loglog(k_axis(2:end), clipv(actual_S4(2:end)), '-', 'LineWidth',1.8, 'Color',[0.47 0.67 0.19]); hold on;
    loglog(k_axis(2:end), clipv(bound1(2:end)), '--', 'LineWidth',1.5, 'Color',[0.85 0.33 0.10]);
    loglog(k_axis(2:end), clipv(bound2(2:end)), '-.', 'LineWidth',1.5, 'Color',[0 0.45 0.74]);
    hold off; grid on;
    xlabel('buoc ngoai k'); ylabel('gia tri (F(x_k)-F^* va cac chan)');
    legend('F(x_k)-F^* thuc te (lich T_k~sqrt(k))', 'chan tu phan tich loai 1', ...
           'chan tu phan tich loai 2', 'Location','southwest', 'FontSize',8);
    title('Kiem chung so hoc cac chan hoi tu');

    filename2 = fullfile(figdir, 'bound-verification.pdf');
    export_fig_pdf(fig2, filename2);

    %% ================= HINH C: suy giam cua eps_k theo cac lich trinh =================
    fig3 = figure('Position',[100 100 640 480]);
    loglog(kk, clipv(R.S1.eps1), '-', 'LineWidth',1.4, 'Color',[0.64 0.08 0.18]); hold on;
    loglog(kk, clipv(R.S2.eps1), '-', 'LineWidth',1.4, 'Color',[0.85 0.33 0.10]);
    loglog(kk, clipv(R.S3.eps1), '-', 'LineWidth',1.4, 'Color',[0.93 0.69 0.13]);
    loglog(kk, clipv(R.S4.eps1), '-', 'LineWidth',1.4, 'Color',[0.47 0.67 0.19]);
    loglog(kk, clipv(R.S5.eps1), '-', 'LineWidth',1.4, 'Color',[0 0.45 0.74]);
    hold off; grid on;
    xlabel('buoc ngoai k'); ylabel('epsilon dang 1 cua k (san tai 10^{-11})');
    legend('T_k=1','T_k=4','T_k ~ log k','T_k ~ sqrt(k)','T_k = T_{ref}', ...
           'Location','southwest', 'FontSize',8);
    title('Do chinh xac loai 1 thuc te dat duoc theo tung lich trinh lap noi');

    filename3 = fullfile(figdir, 'type1-accuracy.pdf');
    export_fig_pdf(fig3, filename3);
end
