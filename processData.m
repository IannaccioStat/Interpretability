function [Lambda_unrot, R_smooth, Z, data_info] = processData(X, Q, opts)
% PROCESSDATA Handles raw data pre-processing, matrix bridge estimation, 
%             PSD spectral projection, normality diagnostics, and unrotated
%             factor extraction (ML or PAF).
%
% Inputs:
%   X    - (N x J) Raw data matrix (observations x variables).
%   Q    - (Scalar) Number of factors to extract.
%   opts - (Struct) Options containing:
%            opts.is_ordinal = 0 (Pearson) or 1 (Polychoric via R bridge)
%
% Outputs:
%   Lambda_unrot - (J x Q) Unrotated factor loading matrix (Lambda_0).
%   R_smooth     - (J x J) Smoothed positive semi-definite correlation matrix.
%   Z            - (N x J) Standardized data matrix.
%   data_info    - (Struct) Combined diagnostic and pipeline metadata.

    if nargin < 3, opts = struct(); end
    
    % Read explicit user flag (Default to 0 / Continuous if not provided)
    is_ordinal_user = get_opt(opts, 'is_ordinal', 0);
    
    [N, J] = size(X);
    data_info = struct();
    data_info.NumVariables = J;
    data_info.NumObservations = N;
    data_info.RequestedFactors = Q;
    
    % Handle constant/zero-variance columns prior to standardizing
    var_x = var(X, 0, 1);
    if any(var_x == 0)
        warning('processData:ZeroVariance', ...
            'Zero-variance variables detected. Standardizing with zero-variance protection.');
        Z = X;
        valid_cols = var_x > 0;
        Z(:, valid_cols) = zscore(X(:, valid_cols));
        Z(:, ~valid_cols) = 0;
    else
        Z = zscore(X);
    end
    
    % ---------------------------------------------------------------------
    % STEP 1: Compute Correlation Matrix R
    % ---------------------------------------------------------------------
    if is_ordinal_user == 0
        % --- PATH 0: Pure Continuous (Pearson) ---
        fprintf('[processData] Mode: Continuous (User Flag = 0). Computing Pearson correlation matrix...\n');
        R = corrcoef(Z);
        
        data_info.CorrelationType  = 'Pearson';
        data_info.NumOrdinalItems  = 0;
        data_info.NumContinuous    = J;
        
    elseif is_ordinal_user == 1
        % --- PATH 1: Ordinal / Likert (Polychoric via R Bridge) ---
        fprintf('[processData] Mode: Ordinal (User Flag = 1).\n');
        
        % Resolve exchange paths
        user_docs = fileparts(userpath); 
        r_folder  = fullfile(user_docs, 'rfiles', 'matlab_bridge');
        
        if ~exist(r_folder, 'dir'), mkdir(r_folder); end
        
        mat_in   = fullfile(r_folder, 'temp_X_data.mat');
        mat_out  = fullfile(r_folder, 'temp_R_poly.mat');
        r_script = fullfile(r_folder, 'calc_polychoric.R');
        
        save(mat_in, 'X', '-v7');
        
        % Locate Rscript.exe dynamically on Windows
        rscript_cmd = 'Rscript';
        if ispc
            r_dirs = dir('C:\Program Files\R\R-*');
            if ~isempty(r_dirs)
                latest_r = r_dirs(end).name;
                rscript_bin = fullfile('C:\Program Files\R', latest_r, 'bin', 'Rscript.exe');
                if exist(rscript_bin, 'file')
                    rscript_cmd = sprintf('"%s"', rscript_bin);
                end
            end
        end
        
        cmd = sprintf('%s "%s"', rscript_cmd, r_script);
        curr_dir = cd(r_folder);
        status = system(cmd);
        cd(curr_dir); % Restore MATLAB path
        
        if status == 0 && exist(mat_out, 'file')
            loaded_data = load(mat_out);
            R = loaded_data.R_poly;
            
            n_ord = get_field_def(loaded_data, 'num_ordinal', J);
            n_con = get_field_def(loaded_data, 'num_continuous', 0);
            
            if isempty(R) || size(R, 1) ~= J || size(R, 2) ~= J
                warning('processData:EmptyRFromBridge', ...
                    'R bridge returned an invalid matrix size. Falling back to Pearson correlation.');
                R = corrcoef(Z);
                n_ord = 0; n_con = J;
            end
            
            data_info.CorrelationType = 'Polychoric (via R psych::polychoric)';
            data_info.NumOrdinalItems = n_ord;
            data_info.NumContinuous   = n_con;
            
            if exist(mat_in, 'file'), delete(mat_in); end
            if exist(mat_out, 'file'), delete(mat_out); end
        else
            error('processData:RBridgeFailed', ...
                'Failed to execute Rscript for Polychoric computation. Verify R installation.');
        end
    end
    
    % ---------------------------------------------------------------------
    % STEP 2: Positive Semi-Definite (PSD) Check & Spectral Projection
    % ---------------------------------------------------------------------
    if any(~isfinite(R(:)))
        warning('processData:NonFiniteValues', 'Non-finite values found in R. Replacing with zeros.');
        R(~isfinite(R)) = 0;
        R(1:J+1:end) = 1;
    end
    
    tol = 1e-6;
    R_sym = (R + R') / 2; % Enforce exact symmetry
    
    [V, D] = eig(R_sym);
    eig_vals = real(diag(D));
    V = real(V);
    
    min_eig = min(eig_vals);
    data_info.MinEigenvalue = min_eig;
    
    if min_eig < tol
        data_info.WasPSDProjected = true;
        fprintf('[processData] Notice: Matrix non-PSD (min eig = %.4e). Projecting to PSD manifold...\n', min_eig);
        
        eig_vals_psd = max(eig_vals, tol);
        R_smooth = V * diag(eig_vals_psd) * V';
        
        % Rescale diagonal to maintain unit variances
        d_scale = sqrt(diag(R_smooth));
        R_smooth = R_smooth ./ (d_scale * d_scale');
        R_smooth(1:J+1:end) = 1; 
    else
        data_info.WasPSDProjected = false;
        R_smooth = R_sym;
        R_smooth(1:J+1:end) = 1;
    end
    
    % ---------------------------------------------------------------------
    % STEP 3: Normality Diagnostics & Factor Extraction
    % ---------------------------------------------------------------------
    fprintf('[processData] Executing normality diagnostics and factor extraction...\n');
    [Lambda_unrot, Psi, diag_info] = select_efa_method(R_smooth, Q, Z, is_ordinal_user);
    
    % Merge diagnostic info into return metadata
    fields = fieldnames(diag_info);
    for k = 1:length(fields)
        data_info.(fields{k}) = diag_info.(fields{k});
    end
    data_info.Uniquenesses = Psi;
end

% =========================================================================
% INTERNAL HELPER 1: SELECT EFA METHOD
% =========================================================================
function [Lambda_unrot, Psi, selection_info] = select_efa_method(R, Q, Z, is_ordinal_user)
    if nargin < 4, is_ordinal_user = 0; end
    [N, J] = size(Z);
    
    % 1. Univariate Skewness and Excess Kurtosis
    sk = skewness(Z, 0);          
    kt = kurtosis(Z, 0) - 3;      
    max_abs_skew = max(abs(sk));
    max_abs_kurt = max(abs(kt));
    
    % 2. Mardia's Normalized Multivariate Kurtosis
    S_cov = cov(Z);
    Z_centered = Z - mean(Z, 1);
    try
        if rcond(S_cov) < 1e-12
            inv_S = pinv(S_cov);
            D_sq = diag(Z_centered * inv_S * Z_centered');
        else
            D_sq = diag(Z_centered * (S_cov \ Z_centered'));
        end
        raw_mardia = sum(D_sq.^2) / N;
        mardia_norm = (raw_mardia - J*(J+2)) / sqrt(8*J*(J+2)/N);
    catch
        mardia_norm = Inf;
    end
    
    % 3. Calculate Normalized Cutoff Bound
    mardia_norm_bound = 0.25 * sqrt(N * J * (J + 2) / 8);
    
    % 4. Decision Rule
    % If explicitly declared ordinal, bypass Gaussian assumptions (non-Gaussian by definition)
    if is_ordinal_user == 1
        is_gaussian = false;
    else
        is_gaussian = (max_abs_skew < 2.0) && ...
                      (max_abs_kurt < 7.0) && ...
                      (abs(mardia_norm) < mardia_norm_bound);
    end
              
    % 5. Execute Factor Extraction
    if is_gaussian
        chosen_method = 'Maximum Likelihood (ML)';
        try
            [Lambda_unrot, Psi] = factoran(Z, Q, 'Rotate', 'none');
            
            if any(Psi <= 0.0051)
                warning('ML hit a boundary constraint (Heywood Case). Falling back to PAF.');
                chosen_method = 'PAF (Fallback from ML Heywood Case)';
                [Lambda_unrot, Psi] = paf(R, Q);
            end
        catch
            warning('ML factoran failed to converge. Falling back to PAF.');
            chosen_method = 'PAF (Fallback from ML Convergence Failure)';
            [Lambda_unrot, Psi] = paf(R, Q);
        end
    else
        chosen_method = 'Principal Axis Factoring (PAF)';
        [Lambda_unrot, Psi] = paf(R, Q);
    end
    
    selection_info = struct();
    selection_info.ChosenMethod             = chosen_method;
    selection_info.IsGaussian               = is_gaussian;
    selection_info.MaxAbsSkewness           = max_abs_skew;
    selection_info.MaxAbsExcessKurtosis     = max_abs_kurt;
    selection_info.MardiaNormalizedKurtosis = mardia_norm;
    selection_info.MardiaNormalizedBound    = mardia_norm_bound;
end

% =========================================================================
% INTERNAL HELPER 2: OPTION PARSERS
% =========================================================================
function val = get_opt(opts, field, default_val)
    if isfield(opts, field) && ~isempty(opts.(field))
        val = opts.(field);
    else
        val = default_val;
    end
end

function val = get_field_def(st, field, default_val)
    if isfield(st, field)
        val = st.(field);
    else
        val = default_val;
    end
end