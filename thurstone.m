function thurstone(Lambda, varargin)
    % THURSTONE Evaluates Thurstone's Simple Structure criteria (1947) 
    % with structured query headers, explicit mathematical descriptions, 
    % and check-level failure tracking.
    
    [J, Q] = size(Lambda);
    
    %% 1. Default Parameters
    default_eps_zero = 0.10;
    default_eps_salient = 0.40;
    
    %% 2. Parse Options (Only essential thresholds configurable)
    p = inputParser;
    addRequired(p, 'Lambda', @isnumeric);
    addParameter(p, 'eps_zero', default_eps_zero, @isnumeric);
    addParameter(p, 'eps_salient', default_eps_salient, @isnumeric);
    
    parse(p, Lambda, varargin{:});
    opts = p.Results;
    
    %% 3. Fixed Thurstone Rules & Base Setup
    % Fixed thresholds derived directly from dimensions
    target_shared = ceil(0.5 * Q);            % Query 4: Majority of required Q negligible loadings
    prop_salient_cross = 1 / (Q + 1);         % Query 5: Inverse factor scaling
    prop_diff = 0.50;                         % Query 3: Absolute majority threshold (50%)
    
    is_zero = abs(Lambda) <= opts.eps_zero;
    is_salient = abs(Lambda) >= opts.eps_salient;
    
    total_salient_vars = sum(any(is_salient, 2));
    max_allowed_cross = ceil(prop_salient_cross * total_salient_vars);
    num_pairs = (Q * (Q - 1)) / 2;
    
    % Per-query trackers for checks evaluated and failures flagged
    checks_per_query = zeros(5, 1);
    fails_per_query  = zeros(5, 1);
    
    fprintf('========================================================================\n');
    fprintf('        THURSTONE SIMPLE STRUCTURE EVALUATION REPORT (J = %d, Q = %d)\n', J, Q);
    fprintf('========================================================================\n\n');
    
    %% =========================================================================
    %% QUERY 1: ROW SIMPLICITY
    %% =========================================================================
    fprintf('------------------------------------------------------------------------\n');
    fprintf('QUERY 1: ROW SIMPLICITY\n');
    fprintf('Description: Each row j (variable) must contain at least one negligible\n');
    fprintf('             loading (|lambda(j,.)| <= %.2f) across all factors.\n', opts.eps_zero);
    fprintf('------------------------------------------------------------------------\n');
    
    checks_per_query(1) = J;
    for j = 1:J
        if sum(is_zero(j, :)) < 1
            fprintf('Row %d FAILED -> no negligible loadings\n', j);
            fails_per_query(1) = fails_per_query(1) + 1;
        end
    end
    if fails_per_query(1) == 0
        fprintf('-> PASSED: All %d rows contain at least one negligible loading.\n', J);
    end
    fprintf('\n');
    
    %% =========================================================================
    %% QUERY 2: COLUMN OVERDETERMINATION
    %% =========================================================================
    fprintf('------------------------------------------------------------------------\n');
    fprintf('QUERY 2: COLUMN OVERDETERMINATION\n');
    fprintf('Description: Each column q (factor) must contain at least Q negligible loadings\n');
    fprintf('             (|lambda(.,q)| <= %.2f) across all variables.\n', opts.eps_zero);
    fprintf('------------------------------------------------------------------------\n');
    
    checks_per_query(2) = Q;
    for q = 1:Q
        zeros_found = sum(is_zero(:, q));
        if zeros_found < Q
            fprintf('Factor %d FAILED -> Not enough negligible loadings (expected: %d, found: %d)\n', ...
                q, Q, zeros_found);
            fails_per_query(2) = fails_per_query(2) + 1;
        end
    end
    if fails_per_query(2) == 0
        fprintf('-> PASSED: All %d factor columns satisfy the overdetermination threshold.\n', Q);
    end
    fprintf('\n');
    
    %% =========================================================================
    %% QUERY 3: FACTOR DIFFERENTIATION (DIRECTIONAL SALIENT SPARSITY)
    %% =========================================================================
    fprintf('------------------------------------------------------------------------\n');
    fprintf('QUERY 3: FACTOR DIFFERENTIATION\n');
    fprintf('Description: For every pair of factors (q,r), a majority (>= 50%%) of the\n');
    fprintf('             variables that are salient (|lambda| >= %.2f) on factor q must\n', opts.eps_salient);
    fprintf('             have negligible loadings (|lambda| <= %.2f) on factor r, AND\n', opts.eps_zero);
    fprintf('             vice versa (r -> q).\n');
    fprintf('------------------------------------------------------------------------\n');
    
    checks_per_query(3) = num_pairs; 
    
    for q = 1:(Q-1)
        for r = (q+1):Q
            pair_failed = false;
            
            % --- Direction 1: q -> r ---
            salient_q_idx = is_salient(:, q);
            salient_q_cnt = sum(salient_q_idx);
            
            if salient_q_cnt > 0
                zeros_r_given_q = sum(salient_q_idx & is_zero(:, r));
                req_zeros_r = ceil(prop_diff * salient_q_cnt);
                
                if zeros_r_given_q < req_zeros_r
                    fprintf('couple (%d,%d) [Dir: %d->%d] FAILED -> Salient on %d: %d; Corresponding zeros on %d: %d (Required majority: %d)\n', ...
                        q, r, q, r, q, salient_q_cnt, r, zeros_r_given_q, req_zeros_r);
                    pair_failed = true;
                end
            end
            
            % --- Direction 2: r -> q ---
            salient_r_idx = is_salient(:, r);
            salient_r_cnt = sum(salient_r_idx);
            
            if salient_r_cnt > 0
                zeros_q_given_r = sum(salient_r_idx & is_zero(:, q));
                req_zeros_q = ceil(prop_diff * salient_r_cnt);
                
                if zeros_q_given_r < req_zeros_q
                    fprintf('couple (%d,%d) [Dir: %d->%d] FAILED -> Salient on %d: %d; Corresponding zeros on %d: %d (Required majority: %d)\n', ...
                        q, r, r, q, r, salient_r_cnt, q, zeros_q_given_r, req_zeros_q);
                    pair_failed = true;
                end
            end
            
            if pair_failed
                fails_per_query(3) = fails_per_query(3) + 1;
            end
        end
    end
    if fails_per_query(3) == 0
        fprintf('-> PASSED: All %d factor pairs meet bidirectional differentiation.\n', checks_per_query(3));
    end
    fprintf('\n');
    
    %% =========================================================================
    %% QUERY 4: SHARED INDEPENDENCE (SHARED ZEROS)
    %% =========================================================================
    fprintf('------------------------------------------------------------------------\n');
    fprintf('QUERY 4: SHARED INDEPENDENCE\n');
    fprintf('Description: For every pair of factors (q,r), the majority of negligible loadings\n');
    fprintf('             (|lambda| <= %.2f) must be shared.\n', opts.eps_zero);
    fprintf('------------------------------------------------------------------------\n');
    
    checks_per_query(4) = num_pairs;
    for q = 1:(Q-1)
        for r = (q+1):Q
            shared_zeros = sum(is_zero(:, q) & is_zero(:, r));
            if shared_zeros < target_shared
                fprintf('couple (%d,%d) FAILED -> expected shared zeros: %d, found: %d\n', ...
                    q, r, target_shared, shared_zeros);
                fails_per_query(4) = fails_per_query(4) + 1;
            end
        end
    end
    if fails_per_query(4) == 0
        fprintf('-> PASSED: All %d factor pairs satisfy the shared zero threshold.\n', num_pairs);
    end
    fprintf('\n');
    
    %% =========================================================================
    %% QUERY 5: LOW COMPLEXITY (MINIMAL CROSS-LOADINGS)
    %% =========================================================================
    fprintf('------------------------------------------------------------------------\n');
    fprintf('QUERY 5: LOW COMPLEXITY\n');
    fprintf('Description: For every pair of factors (q,r), only a minimal number of variables\n');
    fprintf('             may be salient (|lambda| >= %.2f) on both factors at once.\n', opts.eps_salient);
    fprintf('             Threshold: <= %d items (scaled inversely by 1/(Q+1) = %.1f%%\n', ...
        max_allowed_cross, prop_salient_cross * 100);
    fprintf('             across all %d salient variables in the matrix).\n', total_salient_vars);
    fprintf('------------------------------------------------------------------------\n');
    
    checks_per_query(5) = num_pairs;
    for q = 1:(Q-1)
        for r = (q+1):Q
            salient_q = is_salient(:, q);
            salient_r = is_salient(:, r);
            shared_salient = sum(salient_q & salient_r);
            
            if shared_salient > max_allowed_cross
                fprintf('couple (%d,%d) FAILED -> Shared salient loadings (%d/salient on %d, %d/salient on %d) [Found: %d, Max allowed: %d]\n', ...
                    q, r, sum(salient_q), q, sum(salient_r), r, shared_salient, max_allowed_cross);
                fails_per_query(5) = fails_per_query(5) + 1;
            end
        end
    end
    if fails_per_query(5) == 0
        fprintf('-> PASSED: All %d factor pairs satisfy minimal cross-loading limits.\n', num_pairs);
    end
    fprintf('\n');
    
    %% =========================================================================
    %% SUMMARY BREAKDOWN & TOTAL FAILURE PERCENTAGE
    %% =========================================================================
    total_checks = sum(checks_per_query);
    total_fails  = sum(fails_per_query);
    fail_percentage = (total_fails / total_checks) * 100;
    
    fprintf('========================================================================\n');
    fprintf('                       EXPLICIT CHECK BREAKDOWN                         \n');
    fprintf('========================================================================\n');
    fprintf(' Query 1 (Row Simplicity)         : %2d / %2d failed\n', fails_per_query(1), checks_per_query(1));
    fprintf(' Query 2 (Column Overdet.)        : %2d / %2d failed\n', fails_per_query(2), checks_per_query(2));
    fprintf(' Query 3 (Factor Differentiation) : %2d / %2d failed\n', fails_per_query(3), checks_per_query(3));
    fprintf(' Query 4 (Shared Independence)    : %2d / %2d failed\n', fails_per_query(4), checks_per_query(4));
    fprintf(' Query 5 (Low Complexity)         : %2d / %2d failed\n', fails_per_query(5), checks_per_query(5));
    fprintf('------------------------------------------------------------------------\n');
    fprintf(' TOTAL CHECKS EVALUATED          : %d\n', total_checks);
    fprintf(' TOTAL CHECKS FAILED             : %d\n', total_fails);
    fprintf(' OVERALL FAILURE PERCENTAGE      : %.2f%%\n', fail_percentage);
    fprintf('========================================================================\n');
end