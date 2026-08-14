function [Lambda, T] = oblimin(A, gamma, num_starts)
% OBLIMIN Direct oblimin oblique factor rotation using Gradient Projection Algorithm
%         with Multi-Start Random Initialization to avoid local minima.
%
% Syntax:
%   [Lambda, T] = oblimin(A)
%   [Lambda, T] = oblimin(A, gamma)
%   [Lambda, T] = oblimin(A, gamma, num_starts)
%
% Inputs:
%   A          - (p x k) Unrotated factor loadings.
%   gamma      - (Scalar) Parameter: 0 (quartimin), 0.5 (biquartimin), 1 (covarimin). 
%                Default = 0.
%   num_starts - (Optional) Number of random starts. Default = 25.
%
% Outputs:
%   Lambda - (p x k) Rotated factor pattern matrix.
%   T      - (k x k) Transformation matrix compatible with MATLAB conventions 
%            (Lambda = A * T).

    if nargin < 2 || isempty(gamma), gamma = 0; end
    if nargin < 3 || isempty(num_starts), num_starts = 25; end
    
    [p, k] = size(A);
    
    % Kaiser Normalization (row-wise scaling)
    h = sqrt(sum(A.^2, 2));
    h(h < 1e-10) = 1e-10;
    A_norm = A ./ h;
    
    best_f = Inf;
    best_T = eye(k);
    
    % ---------------------------------------------------------------------
    % MULTI-START GRADIENT PROJECTION ALGORITHM (GPA)
    % ---------------------------------------------------------------------
    for s = 1:num_starts
        if s == 1
            T_init = eye(k);
        else
            % Random Orthogonal Matrix via QR Decomposition
            [Q_rand, ~] = qr(randn(k));
            T_init = Q_rand;
        end
        
        [T_opt, f_opt] = run_oblimin_gpa_single(A_norm, T_init, gamma);
        
        % Retain global minimum
        if f_opt < best_f
            best_f = f_opt;
            best_T = T_opt;
        end
    end
    
    T = best_T;
    
    % Compute rotated pattern matrix (un-normalize Kaiser scaling)
    Ti = T \ eye(k);
    Lambda = (A_norm * Ti') .* h;
    
    % Convert T to match MATLAB's "rotatefactors" convention (Lambda = A * T)
    T = (T \ eye(k))'; 
end

% =========================================================================
% INTERNAL HELPER: SINGLE OBLIMIN GPA RUN
% =========================================================================
function [T, f] = run_oblimin_gpa_single(A_norm, T_init, gamma)
    T = T_init;
    [f, G_T] = oblimin_criterion(A_norm, T, gamma);
    
    alpha = 1.0;
    max_iter = 1000;
    ftol = 1e-7;
    
    for iter = 1:max_iter
        % Project gradient onto manifold diag(T' * T) = I
        G_proj = G_T - T * diag(diag(T' * G_T));
        
        if norm(G_proj, 'fro') < ftol, break; end
        
        alpha = alpha * 2;
        for line_step = 1:20
            T_trial = T - alpha * G_proj;
            % Renormalize columns to unit length
            T_trial = T_trial ./ sqrt(sum(T_trial.^2, 1));
            
            [f_trial, G_trial] = oblimin_criterion(A_norm, T_trial, gamma);
            
            if f_trial < f - 0.5 * alpha * sum(G_proj.^2, 'all')
                T = T_trial; 
                f = f_trial; 
                G_T = G_trial;
                break;
            end
            alpha = alpha * 0.5;
        end
    end
end

% =========================================================================
% INTERNAL HELPER: OBLIMIN CRITERION DEFINITION
% =========================================================================
function [f, G_T] = oblimin_criterion(A, T, gamma)
    [p, ~] = size(A);
    
    % Robust inverse computation
    Ti = T \ eye(size(T, 1));
    L = A * Ti';
    L2 = L.^2;
    
    C = eye(p) - (gamma / p) * ones(p, p);
    M = C * L2 * (ones(size(T, 2)) - eye(size(T, 2)));
    
    f = 0.25 * sum(L2 .* M, 'all');
    
    G_L = L .* M;
    G_T = - (L' * G_L * Ti)';
end