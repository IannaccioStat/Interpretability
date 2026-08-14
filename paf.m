function [Lambda, Psi, stats] = paf(Z, Q, max_iter, tol)
% PAF Principal Axis Factoring for Exploratory Factor Analysis.
%
% Syntax:
%   [Lambda, Psi]        = paf(Z, Q)
%   [Lambda, Psi, stats] = paf(Z, Q, max_iter, tol)
%
% Inputs:
%   Z        - (N x J) Data matrix OR (J x J) correlation matrix.
%   Q        - (Scalar) Number of factors to extract.
%   max_iter - (Optional) Maximum iterations (Default: 100).
%   tol      - (Optional) Convergence tolerance for communalities (Default: 1e-5).
%
% Outputs:
%   Lambda   - (J x Q) Unrotated factor loading matrix.
%   Psi      - (J x 1) Unique variances vector.
%   stats    - (Struct) Iteration metadata and convergence status.

    if nargin < 3 || isempty(max_iter), max_iter = 100; end
    if nargin < 4 || isempty(tol),      tol = 1e-5;      end

    % 1. Compute or extract Correlation Matrix
    if size(Z, 1) ~= size(Z, 2)
        R = corr(Z, 'Rows', 'pairwise');
    else
        R = Z;
    end
    
    J = size(R, 1);
    
    % 2. Input Validation
    if Q > J
        error('paf:InvalidQ', 'Number of factors Q (%d) cannot exceed number of variables J (%d).', Q, J);
    end
    
    % Handle Non-Finite Values
    if any(~isfinite(R(:)))
        R(~isfinite(R)) = 0;
        R(1:J+1:end) = 1;
    end
    
    % Force exact matrix symmetry
    R = (R + R') / 2;
    
    % 3. Initial Communality Estimates (SMC with Fallback)
    try
        if rcond(R) < 1e-12
            invR_diag = diag(pinv(R));
        else
            invR_diag = diag(R \ eye(J));
        end
        h2 = 1 - (1 ./ invR_diag);
        
        % Validate SMC estimates
        if any(isnan(h2)) || any(h2 <= 0) || any(h2 >= 1)
            h2 = max(abs(R - eye(J)), [], 2);
        end
    catch
        h2 = max(abs(R - eye(J)), [], 2);
    end
    
    % Clamp initial communalities
    h2 = min(max(h2, 0.01), 0.99);
    
    % 4. Iterative Eigen-decomposition
    R_adj = R;
    converged = false;
    
    for iter = 1:max_iter
        % Replace diagonal with current communalities
        R_adj(1:J+1:end) = h2;
        
        % Spectral Decomposition
        [V, D] = eig(R_adj);
        d = real(diag(D));
        V = real(V);
        
        % Sort eigenvalues in descending order
        [d_sorted, idx] = sort(d, 'descend');
        V_sorted = V(:, idx);
        
        % Extract top Q eigenvalues and eigenvectors
        d_q = max(d_sorted(1:Q), 0); % Truncate negative eigenvalues to 0
        V_q = V_sorted(:, 1:Q);
        
        % Compute factor loadings
        Lambda = V_q * diag(sqrt(d_q));
        
        % Update communalities
        h2_new = sum(Lambda.^2, 2);
        
        % Prevent Heywood cases during iteration
        h2_new = min(max(h2_new, 0.001), 0.999);
        
        % Convergence check
        max_diff = max(abs(h2_new - h2));
        if max_diff < tol
            converged = true;
            break;
        end
        
        h2 = h2_new;
    end
    
    % 5. Calculate Final Uniquenesses (Psi)
    Psi = 1 - sum(Lambda.^2, 2);
    Psi = max(Psi, 0.001); % Enforce positive lower bound
    
    % 6. Populate Diagnostic Stats
    stats = struct();
    stats.Iterations  = iter;
    stats.Converged   = converged;
    stats.MaxAbsDiff  = max_diff;
    stats.FinalH2     = sum(Lambda.^2, 2);
end