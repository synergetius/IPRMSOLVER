function run_iprm_lp()
    lib_path = './netlib/';
    tests = dir(fullfile(lib_path,'*.mat'));
    results = table();
    st = {"QOCO_SOLVED","QOCO_SOLVED_INACCURATE", "QOCO_NUMERICAL_ERROR", "QOCO_MAX_ITER"};
    % test.name = "lp_80bau3b.mat"
    % tests = [test]
    for k = 1:length(tests)
        test = tests(k).name;
        [~, test_name, ~] = fileparts(test);
        load(fullfile(lib_path,test));
        [n, m, p, P, c, A, b, G, h, oscale] = parse_netlib(Problem);

        % n = 1;
        % m = 1;
        % p = 0;
        % P = 0 * speye(1, 1);
        % c = 1;
        % A = sparse([]);
        % b = [];
        % G = -speye(1, 1);
        % h = [-1.25]; 
        solver = iprm;
        settings.verbose = 1;
        settings.max_iters = 200;
        settings.iter_ref_iters = 1;
        settings.epsilon = 1e-12; %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%settings.mu0 = 0.1; %%%%%%%%%%%%%%%%%
        solver.setup(n, m, p, P, c, A, b, G, h, settings);
        out = solver.solve();
        obj = out.obj * oscale
        row = {test_name,st{out.status},out.setup_time_sec,out.solve_time_sec,out.setup_time_sec + out.solve_time_sec, obj, out.phi, out.iters};
        results(k,:) = row;
        results.Properties.VariableNames = ["test","status","setup_time","solve_time","run_time","obj", "phi", "iters"];
        writetable(results, "lp_results.csv");
    end
end
function [n, m, p, P, c, A, b, G, h, oscale] = parse_netlib(Problem)
    A       = Problem.A;
    b       = Problem.b;
    bl      = Problem.aux.lo;
    bu      = Problem.aux.hi;
    c       = Problem.aux.c;
    [p,n]   = size(A);
    scale = 1;
    if scale
        % A * x = b, bl <= x <= bu
        % x乘以cscale，则bl, bu乘以cscale，A 右边乘以Cinv
        % A 左边乘以Rinv，则b除以rscale
        [cscale,rscale] = ScaleData(A,0,0.9);
        C = spdiags(cscale,0,n,n);   Cinv = spdiags(1./cscale,0,n,n);
        R = spdiags(rscale,0,p,p);   Rinv = spdiags(1./rscale,0,p,p);
    
        A  = Rinv*A*Cinv;
        b  = b ./rscale;
        c  = c ./cscale;
        bl = bl.*cscale;
        bu = bu.*cscale;
    end
    fixed   = find(bl==bu);
    blpos   = find(bl< bu & bl>0);
    buneg   = find(bl< bu & bu<0);
    rhs     = b - A(:,fixed)*bl(fixed) ...
        - A(:,blpos)*bl(blpos) ...
        - A(:,buneg)*bu(buneg);

    bscale  = norm(rhs,inf);   
    bscale  = max(bscale,1);
    cscale2  = norm(c,inf);     
    cscale2  = max(cscale2,1);
    oscale = cscale2 * bscale;
    if scale
        % A * x = b, bl <= x <= bu
        % obj = c^T * x
        % x 除以bscale，使b, bl, bu均除以bscale
        % 再让c 除以cscale2，则相应的obj 一共除以了 (cscale2 * bscale)
        % oscale = cscale2 * bscale 这个因子要乘回来
        b       = b /bscale;
        bl      = bl/bscale;
        bu      = bu/bscale;
        c       = c /cscale2;
    end
    G = [speye(n); -speye(n)];
    h = [bu; -bl];
    idx = h ~= inf;
    G = G(idx, :);
    h = h(idx);
    m = size(G, 1);
    P = sparse(n, n);

end