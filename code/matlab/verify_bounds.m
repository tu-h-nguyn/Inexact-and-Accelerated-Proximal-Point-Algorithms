function ok = verify_bounds(seed)
% VERIFY_BOUNDS Kiem chung bang so cac khang dinh ma bao cao chung minh.
%
% Chay khong can giao dien do hoa, in ra bang dat/hong, va tra ma thoat 1 neu
% co bat ky muc nao hong -- de dung truc tiep lam cong kiem tra trong CI.
%
%   octave --no-gui --quiet --eval "exit(~verify_bounds())"
%
% Vi sao kiem tra CAU TRUC chu khong so tung chu so: MATLAB va Octave dung hai
% bo sinh so ngau nhien khac nhau, nen cung mot hat giong van cho hai thuc the
% bai toan khac nhau. So sanh tung chu so se do kim -- no se bao dong ngay khi
% doi moi truong, ma khong he phat hien duoc loi that su. Nhung dieu duoi day
% thi phai dung tren MOI thuc the, vi chung la noi dung cua cac dinh ly.

if nargin < 1 || isempty(seed)
    seed = 7;
end

nfail = 0; ntot = 0;

    function check(label, cond, detail)
        ntot = ntot + 1;
        if cond
            printf('  ok    %s', label);
        else
            printf('  HONG  %s', label);
            nfail = nfail + 1;
        end
        if nargin >= 3 && ~isempty(detail)
            printf('   [%s]', detail);
        end
        printf('\n');
    end

printf('Kiem chung cac chan ly thuyet (seed = %d)\n', seed);

[A, b, ~, par] = setup_problem(seed);
x0 = zeros(par.n, 1);
sched = schedules(par.Tref);

% --- moc tham chieu do chinh xac cao: F*, x* ------------------------------
res_ref = run_outer(A, b, par.rho, par.mu, par.lambda, x0, ...
                    par.K_ref, par.Tref_ref, @(k) par.Tref_ref);
F_star = res_ref.Fvals(end);
x_star = res_ref.xs(:, end);
printf('F* = %.10e\n\n', F_star);

R = struct();
names = fieldnames(sched);
for i = 1:numel(names)
    R.(names{i}) = run_outer(A, b, par.rho, par.mu, par.lambda, x0, ...
                             par.K, par.Tref, sched.(names{i}));
end
Fcl = run_classical(A, b, par.rho, par.mu, par.lambda, x0, par.K, par.Tref);

% =========================================================================
% A. Chan hoi tu cua Dinh ly 3.2 -- khang dinh trung tam
% =========================================================================
printf('A. Chan hoi tu (Dinh ly 3.2): F(x_k)-F* <= beta_k*phi_0 + delta_k\n');
for i = 1:numel(names)
    nm = names{i};
    A0 = R.(nm).Avec(1);
    phi0 = (F_obj(A, b, par.rho, par.mu, x0) - F_star) ...
           + (A0/2) * norm(x_star - x0)^2;
    beta = R.(nm).Avec / A0;
    actual = R.(nm).Fvals - F_star;
    % Dung dung sai tuyet doi nho: F* chi la xap xi so hoc cua gia tri toi uu,
    % nen "vi pham" duoi muc nhieu so hoc khong phai vi pham that.
    tol = 1e-9;
    v1 = sum(actual > beta .* phi0 + R.(nm).delta1 + tol);
    v2 = sum(actual > beta .* phi0 + R.(nm).delta2 + tol);
    check(sprintf('%s: vi pham chan Loai 1', nm), v1 == 0, ...
          sprintf('%d / %d buoc', v1, par.K + 1));
    check(sprintf('%s: vi pham chan Loai 2', nm), v2 == 0, ...
          sprintf('%d / %d buoc', v2, par.K + 1));
end

% =========================================================================
% B. Quan he giua sai so loai 1 va loai 2
% =========================================================================
% CAN THAN: bao cao KHONG khang dinh delta2 <= delta1 tai moi buoc. Cong thuc
% cap nhat khac nhau o cho loai 1 tich luy eta_k = sum eps_i/alpha_i truoc khi
% binh phuong, con loai 2 binh phuong eps_k ngay. Hau qua la loai 2 THIET THOI
% o nhung buoc dau va chi tro nen chat hon ve sau -- do dung la y nghia cua
% "khoi phuc duoc bac hoi tu". Kiem tra dung dieu do, khong kiem tra nhieu hon.
printf('\nB. Quan he giua sai so loai 1 va loai 2\n');
cf2 = (par.lambda*par.rho + 1) / (par.lambda*par.rho);
for i = 1:numel(names)
    nm = names{i};
    check(sprintf('%s: eps2/eps1 = sqrt((lambda*rho+1)/(lambda*rho)) chinh xac', nm), ...
          max(abs(R.(nm).eps2 - R.(nm).eps1 * sqrt(cf2))) < 1e-12);
    % Tai buoc dau tien eta_1 = eps_1/alpha_1 nen (alpha*eta)^2 = eps_1^2, va
    % ty so delta2/delta1 rut gon DUNG BANG (lambda*rho+1)/(lambda*rho).
    check(sprintf('%s: delta2/delta1 tai buoc 1 = %.0f chinh xac', nm, cf2), ...
          abs(R.(nm).delta2(2)/R.(nm).delta1(2) - cf2) < 1e-10, ...
          sprintf('%.10f', R.(nm).delta2(2)/R.(nm).delta1(2)));
    check(sprintf('%s: loai 2 chat hon loai 1 o cuoi day', nm), ...
          R.(nm).delta2(end) < R.(nm).delta1(end), ...
          sprintf('ty so cuoi = %.3g', R.(nm).delta2(end)/R.(nm).delta1(end)));
end
dd = R.S4.delta2 - R.S4.delta1;
cross = find(dd(2:end) <= 0, 1);
check('co dung mot diem giao: loai 2 tu thiet thoi chuyen sang chat hon', ...
      ~isempty(cross) && all(dd(cross+1:end) <= 0), ...
      sprintf('giao tai buoc %d', cross));

% =========================================================================
% C. Tang toc that su xay ra, va no phu thuoc lich trinh lap noi
% =========================================================================
printf('\nC. Tac dung cua tang toc va cua lich trinh lap noi\n');
gap = @(nm) R.(nm).Fvals(end) - F_star;
check('S5 (noi chinh xac) tot hon PPA co dien', ...
      gap('S5') < Fcl(end) - F_star, ...
      sprintf('%.3e  so voi  %.3e', gap('S5'), Fcl(end) - F_star));
check('S1 (T_k=1, cuc tho) pha huy su hoi tu', ...
      gap('S1') > 1e3 * gap('S5'), ...
      sprintf('%.3e  so voi  %.3e', gap('S1'), gap('S5')));
check('lich trinh TANG DAN (S3, S4) giu duoc do chinh xac cua S5', ...
      gap('S3') < 1e3 * gap('S5') && gap('S4') < 1e3 * gap('S5'), ...
      sprintf('S3 %.3e, S4 %.3e, S5 %.3e', gap('S3'), gap('S4'), gap('S5')));
check('lich trinh HANG SO tho (S2) thi khong', ...
      gap('S2') > 1e3 * gap('S5'), ...
      sprintf('S2 %.3e', gap('S2')));
% KHONG do "bac hoi tu" bang do doc hoi quy log-log tren duoi day: voi rho > 0
% bai toan LOI MANH, nen che do tiem can la hoi tu TUYEN TINH (hinh hoc), chu
% khong phai O(1/k^2). Do doc log-log o duoi day chi do cai san sai so so hoc.
% Khang dinh co y nghia la so sanh TRUC TIEP voi hai duong tham chieu o giai
% doan dau, khi bac da tiem can con chua chi phoi.
early = 1:round(par.K/4);
g0 = F_obj(A, b, par.rho, par.mu, x0) - F_star;
check('giai doan dau: S5 nam duoi duong tham chieu O(1/k)', ...
      all(R.S5.Fvals(early+1) - F_star <= g0 ./ early(:)), ...
      sprintf('%d / %d buoc dau', sum(R.S5.Fvals(early+1) - F_star <= g0 ./ early(:)), numel(early)));
check('giai doan dau: PPA co dien cham hon S5 tai moi buoc', ...
      all(Fcl(early+1) >= R.S5.Fvals(early+1) - 1e-12), ...
      sprintf('%d / %d buoc dau', sum(Fcl(early+1) >= R.S5.Fvals(early+1) - 1e-12), numel(early)));

% =========================================================================
% D. Nhat quan noi bo cua vet lap
% =========================================================================
printf('\nD. Nhat quan noi bo\n');
% Gia thiet cua dinh ly: a <= alpha_k^2/((1-alpha_k) A_k lambda_k) <= 2.
% Cong thuc dong kin trong run_outer.m chon alpha la nghiem duong cua
% alpha^2 + (A*lambda)*alpha - (A*lambda) = 0, nen ty so nay bang DUNG 1 --
% nam gon trong cua so cho phep. Kiem tra den chu so may.
rat = R.S4.alpha.^2 ./ ((1 - R.S4.alpha) .* R.S4.Avec(1:end-1) * par.lambda);
check('gia thiet dinh ly: alpha_k^2/((1-alpha_k) A_k lambda) = 1 chinh xac', ...
      max(abs(rat - 1)) < 1e-10, ...
      sprintf('lech lon nhat %.2e', max(abs(rat - 1))));
check('gia thiet dinh ly: ty so nam trong cua so cho phep (0, 2]', ...
      all(rat > 0 & rat <= 2));
check('he thuc truy hoi A_{k+1} = (1-alpha_k) A_k chinh xac', ...
      max(abs(R.S4.Avec(2:end) - (1 - R.S4.alpha) .* R.S4.Avec(1:end-1))) == 0);
check('alpha_k nam trong (0,1) va A_k giam don dieu', ...
      all(R.S4.alpha > 0 & R.S4.alpha < 1) && all(diff(R.S4.Avec) <= 0));
% delta KHONG don dieu -- he so (1-alpha_k) < 1 lam no co the giam. Chi khang
% dinh dieu thuc su dung: khong am.
check('delta1, delta2 khong am tai moi buoc', ...
      all(R.S4.delta1 >= 0) && all(R.S4.delta2 >= 0));
check('S5 dat duoc eps1 nho nhat trong cac lich trinh', ...
      R.S5.eps1(end) <= min([R.S1.eps1(end), R.S2.eps1(end), ...
                             R.S3.eps1(end), R.S4.eps1(end)]) + 1e-12);

printf('\n%d/%d muc dat.\n', ntot - nfail, ntot);
ok = (nfail == 0);
end
