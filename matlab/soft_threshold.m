function x = soft_threshold(v, tau)
% SOFT_THRESHOLD Toan tu Proximal cho chuan L1.
x = sign(v) .* max(abs(v) - tau, 0);
end
