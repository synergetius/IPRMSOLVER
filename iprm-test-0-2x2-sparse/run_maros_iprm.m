function run_maros_iprm()
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
    for k = 1:length(tests)
        test = tests(k).name;
        [~, test_name, ~] = fileparts(test);
        test_name
        mm_data = load("maros/" + test);
        [n, m, p, P, c, A, b, G, h] = parse_maros(mm_data);
        solver = iprm;
        settings.verbose = 1;
        settings.max_iters = 200;
        settings.iter_ref_iters = 1; %%%%%% iter-ref的次数和reg的大小都对求解精度有直接影响
        reg = 1e-8; %%%%%%%%%% 
        settings.kkt_static_reg = reg;
        settings.kkt_dynamic_reg = reg;
        settings.epsilon = 1e-7;
        % settings.epsilon = 1e-5; % 1e-7是QOCO测试采用的精度
        % settings.tau = 0.01; % 影响线搜索的终止条件
        % settings.sigma = 0.01; % 影响rho的更新
        % settings.delta = 0.5; % 线搜索的倍率
        % settings.gamma0 = 0.97;%0.97; % 影响很明显。论文使用gamma0 = 0.001
        % settings.max_iters = 200; % 10000;
        solver.setup(n, m, p, P, c, A, b, G, h, settings);
        out = solver.solve();
        out.phi
        row = {test_name,st{out.status},out.setup_time_sec,out.solve_time_sec,out.setup_time_sec + out.solve_time_sec, out.obj, out.phi, out.iters};
        results(k,:) = row;
        results.Properties.VariableNames = ["test","status","setup_time","solve_time","run_time","obj", "phi", "iters"];
        writetable(results, "maros_results.csv");
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

