function [pass] = no_constraints()
    n = 6;
    m = 0;
    p = 0;

    P = diag([1,2,3,4,5,6]);
    c = [1;2;3;4;5;6];
    A = [];
    b = [];
    G = [];
    h = [];

    solver = iprm;
    solver.setup(n, m, p, P, c, A, b, G, h)
    out = solver.solve();

    tol = 1e-3;
    x_exp = [-1, -1, -1, -1, -1, -1];
    %out.status, norm(out.x - x_exp) 
    %pass = norm(out.x - x_exp) < tol;
    pass = (out.status == 1) && norm(out.x - x_exp) < tol;
end
