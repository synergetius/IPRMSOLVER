function run_dimacs_debug()
    results = table();

    st = {"QOCO_SOLVED","QOCO_SOLVED_INACCURATE", "QOCO_NUMERICAL_ERROR", "QOCO_MAX_ITER"};
    test.name = "nb.mat" % 是最小的一个数据集，无缩放时数值正常但解不出来
    %test.name = "qssp30.mat" % 可以解出
    tests = [test]
    for k = 1:length(tests)
        test = tests(k).name;
        [~, test_name, ~] = fileparts(test);
        test_name
        mm_data = load("dimacs-sedumi/" + test);
        [n, m, p, P, c, A, b, G, h, l, nsoc, q, oscale] = parse_dimacs(mm_data);
        disp({n, m, p, P, c, A, b, G, h, l, nsoc, q, oscale})
        %%%%%%%%% 要确定传入的数组形状合适
    
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
        writetable(results, "dimacs_results.csv");
    end
end
function [n, m, p, P, c, A, b, G, h, l, nsoc, q, oscale] = parse_dimacs(mm_data)
    if isfield(mm_data, 'At')
        mm_data.A = mm_data.At';
    end
    A = mm_data.A;
    b = full(mm_data.b);
    %b
    c = full(mm_data.c);
    K = mm_data.K;
    [p, n] = size(A);
    G = -speye(n);
    h = zeros(n, 1);
    m = n;
    P = sparse(n, n);
    l = K.l;
    q = full(K.q)';
    nsoc = length(q);
    oscale = 1; %% 无缩放
end