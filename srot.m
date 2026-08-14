function [Lambda, T, f_opt] = srot(A, method_name, num_starts)
% SROT Oblique factor rotation using sparse simplicity criteria with 
%      Multi-Start Gradient Projection Algorithm (GPA) to avoid local minima.
%
% Syntax:
%   [Lambda, T] = srot(A, method_name)
%   [Lambda, T] = srot(A, method_name, num_starts)
%   [Lambda, T, f_opt] = srot(...)
%
% Inputs:
%   A           - (p x k) Unrotated factor loading matrix (p variables, k factors).
%   method_name - (String) Rotation criterion:
%                 'GEOMIN'     - Geomin criterion (Yates, 1987).
%                 'INFOMAX'    - Infomax criterion (McKeon, 1968).
%                 'MCCAMMON'   - Minimum Entropy criterion (McCammon, 1966).
%                 'CLF_SPARSE' - Component Loss Function / L0-approximation 
%                                (Jennrich, 2006).
%   num_starts  - (Optional) Number of random initializations. Default: 25.
%
% Outputs:
%   Lambda - (p x k) Rotated factor pattern loading matrix.
%   T      - (k x k) Direct transformation matrix such that Lambda = A_norm * T_inv_t.
%   f_opt  - (Scalar) Final value of the simplicity objective function.
%
% References:
%   - Yates, A. (1987). Multivariate Exploratory Data Analysis. SUNY Press.
%   - McKeon, J. J. (1968). Infomax: A rotation criterion for max entropy. 
%     Psychometrika, 33(2), 219-228.
%   - McCammon, R. B. (1966). Principal component analysis and its application 
%     in large-scale correlation studies. Journal of Geology, 74(5), 721-733.
%   - Jennrich, R. I. (2006). Rotation to simple structure using component loss 
%     functions. Psychometrika, 71(1), 173-191.
%   - Bernaards, C. A., & Jennrich, R. I. (2005). Gradient projection algorithms 
%     for rotation in factor analysis. Educational and Psychological Measurement, 
%     65(5), 766-796.

    %% 1. Input Validation
    narginchk(2, 3);
    validateattributes(A, {'numeric'}, {'2d', 'real', 'finite', 'nonempty'}, mfilename, 'A', 1);
    
    [p, k] = size(A);
    if k < 2
        warning('srot:SingleFactor', 'Matrix A has less than 2 factors. Rotation skipped.');
        Lambda = A;
        T = 1;
        f_opt = 0;
        return;
    end
    
    if nargin < 3 || isempty(num_starts)
        num_starts = 25;
    end
    
    valid_methods = {'GEOMIN', 'INFOMAX', 'MCCAMMON', 'CLF_SPARSE'};
    method_name = upper(string(method_name));
    if ~ismember(method_name, valid_methods)
        error('srot:InvalidMethod', ...
            'Unknown method ''%s''. Supported methods: %s', method_name, strjoin(valid_methods, ', '));
    end

    %% 2. Kaiser Normalization
    h2 = sum(A.^2, 2);
    h = sqrt(h2);
    h(h < 1e-10) = 1e-10;
    A_norm = A ./ h;

    %% 3. Multi-Start Gradient Projection Optimization
    best_f = Inf;
    best_T = eye(k);

    for s = 1:num_starts
        if s == 1
            T_init = eye(k);
        else
            % Random orthogonal matrix initialization via QR decomposition
            [Q_rand, ~] = qr(randn(k));
            T_init = Q_rand;
        end

        [T_opt, f_opt_run] = run_gpa_single(A_norm, T_init, method_name);

        if f_opt_run < best_f
            best_f = f_opt_run;
            best_T = T_opt;
        end
    end

    %% 4. Final Solution Assembly & Standard Normalization
    % Normalize columns of T so factor correlation matrix Phi has unit diagonal
    col_norms = sqrt(sum(best_T.^2, 1));
    best_T = best_T ./ col_norms;
    
    % Compute primary pattern loading matrix: Lambda = A_norm * inv(T)'
    Ti = best_T \ eye(k);
    Lambda_norm = A_norm * Ti';
    
    % Un-normalize Kaiser scaling
    Lambda = Lambda_norm .* h;
    T = best_T;
    f_opt = best_f;
end

% =========================================================================
% INTERNAL HELPER: SINGLE GPA OPTIMIZATION RUN
% =========================================================================
function [T, f] = run_gpa_single(A_norm, T_init, method_name)
    T = T_init;
    
    % Safe initial evaluation
    [f, G_T] = sparse_criterion(A_norm, T, method_name);
    if isinf(f) || isnan(f)
        f = Inf;
        return;
    end

    alpha = 1.0;
    max_iter = 1000;
    ftol = 1e-7;

    for iter = 1:max_iter
        % Project gradient onto manifold of unit column-length matrices
        G_proj = G_T - T * diag(sum(T .* G_T, 1));
        
        proj_norm = sqrt(sum(G_proj.^2, 'all'));
        if proj_norm < ftol
            break; 
        end

        % Line Search (Armijo condition)
        alpha = alpha * 2;
        stepped = false;
        
        for line_step = 1:20
            T_trial = T - alpha * G_proj;
            col_norms = sqrt(sum(T_trial.^2, 1));
            col_norms(col_norms < eps) = eps;
            T_trial = T_trial ./ col_norms;
            
            [f_trial, G_trial] = sparse_criterion(A_norm, T_trial, method_name);

            if f_trial < f - 0.5 * alpha * sum(G_proj.^2, 'all')
                T = T_trial; 
                f = f_trial; 
                G_T = G_trial;
                stepped = true;
                break;
            end
            alpha = alpha * 0.5;
        end
        
        if ~stepped
            break; % Convergence stalling guard
        end
    end
end

% =========================================================================
% INTERNAL HELPER: SPARSE CRITERIA & ANALYTICAL GRADIENTS
% =========================================================================
function [f, G_T] = sparse_criterion(A, T, method)
    [p, k] = size(A);
    
    % Compute L = A * inv(T)' without explicit matrix inversion
    if rcond(T) < 1e-12
        f = Inf;
        G_T = zeros(size(T));
        return;
    end
    
    Ti = T \ eye(k);
    L = A * Ti';
    L2 = L.^2;
    eps_val = 0.01;

    switch method
        case 'GEOMIN'
            % Geomin criterion (Yates, 1987)
            geom = prod(L2 + eps_val, 2).^(1/k);
            f = sum(geom);
            G_L = (2/k) * L .* (geom ./ (L2 + eps_val));

        case 'INFOMAX'
            % Infomax entropy minimization (McKeon, 1968)
            S_L = max(sum(L2, 'all'), eps);
            P = L2 / S_L;
            P = max(P, 1e-12);
            f = sum(P .* log(P), 'all');
            G_L = (2 / S_L) * L .* (1 + log(P) - f);

        case 'MCCAMMON'
            % McCammon Minimum Entropy (McCammon, 1966)
            col_sums = max(sum(L2, 1), eps);
            P = L2 ./ col_sums;
            P = max(P, 1e-12);
            
            entropy_cols = sum(P .* log(P), 1);
            f = -sum(entropy_cols);
            
            G_L = -2 * (L ./ col_sums) .* (1 + log(P) - entropy_cols);

        case 'CLF_SPARSE'
            % Concave Component Loss / Hard Sparsity L0-approx (Jennrich, 2006)
            exp_term = exp(-L2 / (2 * eps_val^2));
            f = sum(1 - exp_term, 'all');
            G_L = (L / (eps_val^2)) .* exp_term;
    end

    % Chain rule: Transform gradient w.r.t loadings (G_L) to gradient w.r.t T (G_T)
    G_T = -(L' * G_L * Ti)';
end