function run_scale_maros()
    tests = dir(fullfile('maros/','*.mat'));
    results = table();

    st = {"QOCO_SOLVED","QOCO_SOLVED_INACCURATE", "QOCO_NUMERICAL_ERROR", "QOCO_MAX_ITER"};
    for k = 1:length(tests)
        test = tests(k).name;
        [~, test_name, ~] = fileparts(test);
        test_name
        mm_data = load("maros/" + test);
        [n, m, p, P, c, A, b, G, h, l, nsoc, q, oscale] = parse_maros(mm_data);
        solver = iprm;
        settings.verbose = 0;
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
        writetable(results, "scale_maros_results.csv");
    end
end
function [n, m, p, P, c, Aeq, beq, G, h, l, nsoc, q, oscale] = parse_maros(mm_data)
    % 提取变量数量
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
    
    % 数据缩放处理

    % 对约束矩阵 Amm 和右侧向量 rl 进行缩放
    [cscale, rscale] = ScaleData(Amm, 0, 0.9);
    C = spdiags(cscale, 0, n, n);
    Cinv = spdiags(1 ./ cscale, 0, n, n);
    R = spdiags(rscale, 0, size(Amm, 1), size(Amm, 1));
    Rinv = spdiags(1 ./ rscale, 0, size(Amm, 1), size(Amm, 1));
    
    % rl <= Amm x <= ru, lb <= x <= ub
    % x 乘以cscale 则lb, ub 乘以cscale, c除以cscale, P 左右乘上Cinv' 和Cinv，Amm 右边乘以Cinv
    % 对rl <= Amm x <= ru整体再除以rscale，则Amm左边乘以Rinv，rl, ru除以rscale
    Amm = Rinv * Amm * Cinv;
    rl = rl ./ rscale;
    ru = ru ./ rscale;
    lb = lb .* cscale;
    ub = ub .* cscale;
    c = c ./ cscale;
    P = Cinv' * P * Cinv;


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
    
    %
    fixed = find(lb == ub);
    blpos = find(lb < ub & lb > 0);
    buneg = find(lb < ub & ub < 0);

    rhs = beq - Aeq(:, fixed) * lb(fixed) ...
        - Aeq(:, blpos) * lb(blpos) ...
        - Aeq(:, buneg) * ub(buneg);
    % 
    bscale = norm(rhs, inf);
    bscale = max(bscale, 1);
    cscale2 = norm(c, inf);
    cscale2 = max(cscale2, 1);
    %
    % 进行缩放：
    % (1) 将x除以bscale，则Aeq x = beq中beq 要除以bscale，lb <= x <= ub 中lb,
    % ub要除以bscale，lineq <= Aineq x <= uineq 中 lineq, uineq 要除以bscale
    % c 要乘以bscale, P 要乘以bscale^2
    % (2) 对目标函数中c^T * x 的项中的c 再除以(bscale * cscale2) （这样总体上c 变为c / cscale2），
    % 为了维持目标函数一致，则P应该再除以 (bscale * cscale2)，则总体上P变为 P * bscale / cscale2
    beq = beq / bscale;
    lb = lb / bscale;
    ub = ub / bscale;
    lineq = lineq / bscale;
    uineq = uineq / bscale;
    c = c / cscale2;
    P = P * bscale / cscale2;
    % 总体上，目标函数被除以了bscale * cscale2，它就是oscale（objective function的缩放值）
    oscale = bscale * cscale2;

    % 重新构建不等式约束矩阵 G 和向量 h
    G = [speye(n); -speye(n); Aineq; -Aineq];
    h = [ub; -lb; uineq; -lineq];
    
    % 移除无穷值
    idx = h ~= inf;
    G = G(idx, :);
    h = h(idx);
    
    % 计算不等式约束的数量和等式约束的数量
    m = size(G, 1);
    p = size(Aeq, 1);
   
    l = m;
    nsoc = 0;
    q = [];
   
end