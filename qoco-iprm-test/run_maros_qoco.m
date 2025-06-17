function run_maros_qoco()
    tests = dir(fullfile('maros/','*.mat'));
    results = table();
    
    %results = table([],[],[],[],[],[],[], "VariableNames",{"test","status","setup_time","solve_time","run_time","obj","iters"});
    
    %test = {};
    %status = {};
    %setup_time = {};
    %solve_time = {};
    %run_time = {};
    %obj = {};
    %iters = {};
    %results = table(test, status, setup_time, solve_time, run_time, obj, iters);
    %colNames = {"test", "status", "setup_time", "solve_time", "run_time", "obj", "iters"};
    %colTypes = {'string', 'string', 'double', 'double', 'double', 'double', 'double'};
    %results = preallocate('table', colNames, colTypes);
    st = {"QOCO_SOLVED","QOCO_SOLVED_INACCURATE", "QOCO_NUMERICAL_ERROR", "QOCO_MAX_ITER"};
    test.name = "QCAPRI.mat";
    tests = [test];
    for k = 1:length(tests)
        test = tests(k).name;
        [~, test_name, ~] = fileparts(test);
        test_name
        mm_data = load("maros/" + test);
        [n, m, p, P, c, A, b, G, h] = parse_maros(mm_data);
        l = m;
        nsoc = 0;
        q = [];
        solver = qoco;
        settings.verbose = 1;
        %settings.epsilon = 0.1;
        settings.max_iters = 200;
        solver.setup(n, m, p, P, c, A, b, G, h, l, nsoc, q, settings);
        out = solver.solve();
        %out.phi
        %row = {test_name,st{out.status},out.setup_time_sec,out.solve_time_sec,out.setup_time_sec + out.solve_time_sec, out.obj, out.phi, out.iters};
        %results(k,:) = row;
        %results.Properties.VariableNames = ["test","status","setup_time","solve_time","run_time","obj", "phi", "iters"];
        %writetable(results, "maros_results.csv");
    end
    %mm_data = load("maros/AUG2D.mat");
    
    % 考虑要记录哪些测试信息
end
function [n, m, p, P, c, Aeq, beq, G, h] = parse_maros(mm_data)
    n = length(mm_data.xl);
    % 提取二次项系数矩阵
    P = mm_data.H;
    
    % 提取线性项系数向量
    c = squeeze(mm_data.g);
    
    % 提取约束矩阵和向量
    Amm = mm_data.A;
    rl = squeeze(mm_data.al);
    ru = squeeze(mm_data.au);
    lb = squeeze(mm_data.xl);
    ub = squeeze(mm_data.xu);
    
    % 区分等式约束和不等式约束
    eq_idx = find(rl == ru);
    ineq_idx = find(rl ~= ru);
    
    % 提取等式约束
    Aeq = Amm(eq_idx, :);
    beq = rl(eq_idx);
    
    % 提取不等式约束
    Aineq = Amm(ineq_idx, :);
    uineq = ru(ineq_idx);
    lineq = rl(ineq_idx);
    
    % 构建不等式约束矩阵 G 和向量 h
    G = [speye(n); -speye(n); Aineq; -Aineq];
    h = [ub; -lb; uineq; -lineq];
    
    % 移除无穷值
    idx = h ~= inf;
    G = G(idx, :);
    h = h(idx);
    
    % 计算不等式约束的数量和等式约束的数量
    m = size(G, 1);
    p = size(Aeq, 1);
    
    
end

