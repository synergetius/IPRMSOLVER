function run_dimacs()
    tests = dir(fullfile('dimacs-sedumi/','*.mat'));
    results = table();

    st = {"QOCO_SOLVED","QOCO_SOLVED_INACCURATE", "QOCO_NUMERICAL_ERROR", "QOCO_MAX_ITER"};
    %test.name = "nql30new.mat" % 是最小的一个数据集，无缩放时数值正常但解不出来
    %test.name = "qssp30.mat" % 可以解出
    % sched_100_100_orig
    %tests = [test]
    for k = 1:length(tests)
        test = tests(k).name;
        [~, test_name, ~] = fileparts(test);
        test_name
        mm_data = load("dimacs-sedumi/" + test);
        [n, m, p, P, c, A, b, G, h, l, nsoc, q, oscale] = parse_dimacs(mm_data);
        %%%%%%%%% 要确定传入的数组形状合适
        
        disp(size(P))
        disp(size(b))
        disp(size(c))
        disp(size(q))
        solver = iprm;
        settings.verbose = 1;
        settings.max_iters = 200;
        settings.iter_ref_iters = 1;
        %%%%%
        %settings.linesearch_iters = 10; %%%%%%%%%%%%%%%% 测试
        %settings.delta = 0.4;

        %settings.linesearch_iters = 7; %%%%%%%%%%%%%%%% 测试
        %settings.delta = 0.3;
        %%%%%%%
        settings.epsilon = 1e-8;
        reg = 1e-8;
        settings.kkt_static_reg = reg;
        settings.kkt_dynamic_reg = reg;
        % settings.epsilon = 1e-8;
        % settings.rho0 = 1e-3;
        % settings.delta = 0.5;
        % settings.sigma = 0.01;
        % settings.tau = 0.01;
        % settings.mu0 = 0.1;
        % settings.eta = 10;
        % settings.gamma0 = 0.1;
        %%%%%%%%
        solver.setup(n, m, p, P, c, A, b, G, h, l, nsoc, q, settings);
        out = solver.solve();
        obj = out.obj * oscale;
        out.phi, obj
        row = {test_name,st{out.status},out.setup_time_sec,out.solve_time_sec,out.setup_time_sec + out.solve_time_sec, obj, out.phi, out.iters};
        if out.status ~= 1
            row{1, 5} = inf;
            row{1, 8} = inf;
        end
        results(k,:) = row;
        
        results.Properties.VariableNames = ["test","status","setup_time","solve_time","run_time","obj", "phi", "iters"];
        writetable(results, "scale_dimacs_results.csv");
    end
end
function [n, m, p, P, c, A, b, G, h, l, nsoc, q, oscale] = parse_dimacs(mm_data)
    if isfield(mm_data, 'At')
        mm_data.A = mm_data.At';
    end
    A = mm_data.A;
    b = full(mm_data.b);
    c = full(mm_data.c);
    K = mm_data.K;
    [p, n] = size(A);
    
    scale = 1;
    if scale
        % A * x = b, x \in K
        % x乘以cscale，则A 右边乘以Cinv
        % A 左边乘以Rinv，则b除以rscale
        [cscale,rscale] = ScaleData(A,0,0.9);
        C = spdiags(cscale,0,n,n);   Cinv = spdiags(1./cscale,0,n,n);
        R = spdiags(rscale,0,p,p);   Rinv = spdiags(1./rscale,0,p,p);
    
        A = Rinv*A*Cinv;
        b = b ./rscale;
        c = c ./cscale;

        rhs = b;
        bscale  = norm(rhs,inf);   
        bscale  = max(bscale,1);
        cscale2  = norm(c,inf);     
        cscale2  = max(cscale2,1);
        oscale = cscale2 * bscale;
        % A * x = b, x \in K
        % obj = c^T * x
        % x 除以bscale，使b除以bscale
        % 再让c 除以cscale2，则相应的obj 一共除以了 (cscale2 * bscale)
        % oscale = cscale2 * bscale 这个因子要乘回来
        b = b /bscale;
        c = c /cscale2;
    end
    G = -speye(n);
    h = zeros(n, 1);
    m = n;
    P = sparse(n, n);
    l = K.l;
    q = full(K.q)';
    nsoc = length(q);

end