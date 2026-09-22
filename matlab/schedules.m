function sched = schedules(Tref)
% SCHEDULES Cac lich trinh so vong lap noi T_k.
%
% Day la bien doc lap cua ca bai toan: thuat toan diem gan ke khong chinh xac
% cho phep giai bai toan con mot cach xap xi, va cau hoi la xap xi TOI DAU thi
% van giu duoc bac hoi tu O(1/k^2).
sched = struct( ...
    'S1', @(k) 1, ...                          % Hang so, cuc tho
    'S2', @(k) 4, ...                          % Hang so, tho
    'S3', @(k) round(2*log2(k+2)) + 1, ...     % Tang cham (logarit)
    'S4', @(k) round(1.5*sqrt(k+1)) + 2, ...   % Tang vua (can bac hai)
    'S5', @(k) Tref);                          % Chinh xac (tham chieu)
end
