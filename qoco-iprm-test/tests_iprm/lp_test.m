function [pass] = lp_test()
    n = 2;
    m = 3;
    p = 0;
    P = [];
    c = [-1;-2];
    A = [];
    b = [];
    G = [-1 0 0;0 -1 0;1 1 0];
    h = [0;0;1];
    solver = iprm;
    solver.setup(n, m, p, P, c, A, b, G, h)
    out = solver.solve();

    tol = 1e-3;
    x_exp = [0.0;1.0];
    xi_exp = [0.0;1.0;0.0];
    s_exp = [1.0;0.0;2.0];
    out.status, norm(out.x - x_exp), norm(out.xi - xi_exp), norm(out.s - s_exp)
    %%%% status 还没设置好 （成功收敛status = 1）
    pass = (out.status == 1) &&  norm(out.x - x_exp) < tol && ...
        norm(out.xi - xi_exp) < tol && norm(out.s - s_exp) < tol;
end
