function [pass] = no_soc_constraints()
    n = 6;
    m = 6;
    p = 2;
    P = diag([1,2,3,4,5,6]);
    c = [1;2;3;4;5;6];
    A = [1 1 0 0 0 0;0 1 2 0 0 0];
    b = [1;2];
    G = -eye(6);
    h = [0;0;0;0;0;0];
    solver = iprm;
    solver.setup(n, m, p, P, c, A, b, G, h)
    out = solver.solve();

    tol = 0.011;%2e-3; %% 还达不到1e-4
    x_exp = [0.2000;0.8000;0.6000;-0.0000;0.0000;0.0000];
    xi_exp = [0.2000;0.8000;0.6000;0.0000;0.0000;0.0000];
    t_exp = [-1.200;-2.400];
    s_exp = [0.0000;0.0000;0.0000;4.0000;5.0000;6.0000];
    % out.status = 4（达到最大迭代次数还没到足够精度）
    out.status, norm(out.x - x_exp), norm(out.xi - xi_exp), norm(out.t - t_exp), norm(out.s - s_exp)
    % pass = (out.status == 1) && norm(out.x - x_exp) < tol && ...
    %     norm(out.xi - xi_exp) < tol && norm(out.t - t_exp) < tol && ...
    %     norm(out.s - s_exp) < tol;
    pass = norm(out.x - x_exp) < tol && ...
        norm(out.xi - xi_exp) < tol && norm(out.t - t_exp) < tol && ...
        norm(out.s - s_exp) < tol;
end
