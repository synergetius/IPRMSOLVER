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
    

    settings.verbose = 1;
    settings.epsilon = 1e-7;
    settings.iter_ref_iters = 1; 
    settings.tau = 1e-2; % 影响线搜索的终止条件
    settings.sigma = 0.01; % 影响rho的更新
    settings.delta = 0.5; % 线搜索的倍率 %%%%%%%%%%%%%%%%%%%%%
    settings.gamma0 = 0.1; % 影响很明显
    settings.max_iters = 200; % 10000;
    % 问题现象描述：gap = z - xi 维持在1e-5无法下降到1e-7，从几次迭代开始，每次迭代也找不到合适的步长
    % 找不到合适的步长：要么是步长的计算出了问题，要么是迭代方向的计算出了问题
    % 发现新的现象：在gap达到4e-5后，Dxi一下子降到1e-7以下的数量级，这直接导致xi的迭代和gap的下降会变得非常缓慢。Dxi应该已经出了问题
    solver.setup(n, m, p, P, c, A, b, G, h, settings);
    out = solver.solve();

    tol = 1e-4; %% 还达不到1e-4
    x_exp = [0.2000;0.8000;0.6000;-0.0000;0.0000;0.0000];
    xi_exp = [0.2000;0.8000;0.6000;0.0000;0.0000;0.0000];
    t_exp = [-1.200;-2.400];
    s_exp = [0.0000;0.0000;0.0000;4.0000;5.0000;6.0000];
    % out.status = 4（达到最大迭代次数还没到足够精度）
    out.status, out.phi, out.iters
    % "x", out.x
    % "xi", out.xi
    % "t", out.t
    % "s", out.s

    %out.status, norm(out.x - x_exp), norm(out.xi - xi_exp), norm(out.t - t_exp), norm(out.s - s_exp)
    pass = (out.status == 1) && norm(out.x - x_exp) < tol && ...
        norm(out.xi - xi_exp) < tol && norm(out.t - t_exp) < tol && ...
        norm(out.s - s_exp) < tol;
end
