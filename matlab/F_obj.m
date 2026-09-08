function v = F_obj(A, b, rho, mu, x)
% F_OBJ Tinh gia tri ham muc tieu F(x).
v = 0.5*norm(A*x - b)^2 + (rho/2)*norm(x)^2 + mu*norm(x, 1);
end
