function report = interp(Lambda, Phi, opts)
% INTERP Quantitative assessment of factor model interpretability.
%
% Syntax:
%   report = interp(Lambda)
%   report = interp(Lambda, Phi)
%   report = interp(Lambda, Phi, opts)
%
% Description:
%   Evaluates the overall interpretability and structural properties of a 
%   factor analysis solution based on three core dimensions:
%     1. Explained Common Variance (ECV) - Variance captured by factors.
%     2. Generalized Column Sparseness (GS_w) - Norm-ratio metric measuring 
%        simple structure (sparsity) per factor.
%     3. Orthogonality Index - Root-Mean-Square (RMS) distance of factor 
%        correlation matrix Phi from an identity matrix.
%
% Inputs:
%   Lambda - (J x Q) Factor pattern loading matrix.
%   Phi    - (Q x Q) Factor correlation matrix. Optional (Default: eye(Q)).
%   opts   - (Struct) Optional parameters to configure evaluation:
%            .p1     - Numerator norm degree for sparseness (Default: 2).
%            .p2     - Denominator norm degree for sparseness (Default: 4).
%            .weights- Weighting scheme: 'orthogonalized' (Default), 
%                      'structure', or 'direct'.
%            .Method - String label for extraction method (e.g., 'ML').
%            .Rotation- String label for rotation method (e.g., 'Varimax').
%
% Outputs:
%   report - (Struct) Contains full metrics, factor details, and diagnostics.

    %% 1. Input Parsing and Default Setup
    narginchk(1, 3);
    
    validateattributes(Lambda, {'numeric'}, {'2d', 'nonempty', 'real', 'finite'}, mfilename, 'Lambda', 1);
    [J, Q] = size(Lambda);
    
    if nargin < 2 || isempty(Phi)
        Phi = eye(Q);
    else
        validateattributes(Phi, {'numeric'}, {'square', 'real', 'finite'}, mfilename, 'Phi', 2);
        if size(Phi, 1) ~= Q
            error('interp:DimensionMismatch', ...
                'Dimension mismatch: Phi must be (%d x %d) to match Lambda columns.', Q, Q);
        end
    end
    
    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    
    p1            = get_opt(opts, 'p1', 2);
    p2            = get_opt(opts, 'p2', 4);
    weight_scheme = lower(get_opt(opts, 'weights', 'orthogonalized'));
    method_label  = get_opt(opts, 'Method', 'Unspecified');
    rot_label     = get_opt(opts, 'Rotation', 'Unspecified');
    
    %% 2. Positive Semi-Definite (PSD) Check & Symmetric Matrix Square Root
    tol = 1e-6;
    Phi_sym = (Phi + Phi') / 2; % Ensure exact symmetry
    [V_phi, D_phi] = eig(Phi_sym);
    eig_vals = diag(D_phi);
    
    was_projected = false;
    if any(eig_vals < tol)
        was_projected = true;
        warning('interp:NonPSD', ...
            'Factor correlation matrix Phi is not positive semi-definite (min eigenvalue = %.4e). Projecting to nearest PSD matrix.', ...
            min(eig_vals));
        
        eig_vals_psd = max(eig_vals, tol);
        Phi_sym = V_phi * diag(eig_vals_psd) * V_phi';
        
        d_scale = sqrt(diag(Phi_sym));
        Phi_sym = Phi_sym ./ (d_scale * d_scale');
        
        % Re-eigen-decompose post-projection
        [V_phi, D_phi] = eig(Phi_sym);
        eig_vals = max(diag(D_phi), 0);
    end
    
    Phi = Phi_sym;
    Phi_sqrt = V_phi * diag(sqrt(max(eig_vals, 0))) * V_phi';
    
    %% 3. Calculate Explained Common Variance (ECV)
    Sigma_hat_common = Lambda * Phi * Lambda';
    ECV = trace(Sigma_hat_common) / J;
    
    %% 4. Vectorized Column-wise Generalized Sparseness (gs_vec)
    norm_const = J^((1/p1) - (1/p2)); % Pre-factor for Lp1/Lp2 generalization
    
    norm_p1 = vecnorm(Lambda, p1, 1)';
    norm_p2 = vecnorm(Lambda, p2, 1)';
    
    gs_vec = zeros(Q, 1);
    valid_mask = (norm_p2 > eps);
    
    if norm_const > 1
        ratios = norm_p1(valid_mask) ./ norm_p2(valid_mask);
        gs_vec(valid_mask) = (norm_const - ratios) / (norm_const - 1);
    end
    
    %% 5. Compute Weight Schemes
    [w_str, w_orth, w_dir] = compute_all_weights(Lambda, Phi, Phi_sqrt);
    
    switch weight_scheme
        case 'orthogonalized'
            w_active = w_orth;
        case 'structure'
            w_active = w_str;
        case 'direct'
            w_active = w_dir;
        otherwise
            error('interp:InvalidScheme', ...
                'Invalid weighting scheme ''%s''. Choose ''orthogonalized'', ''structure'', or ''direct''.', weight_scheme);
    end
    
    GS_w = sum(w_active .* gs_vec);
    I_index = sqrt(ECV * GS_w);
    
    %% 6. Orthogonality Check (RMS Distance from Identity Matrix)
    if Q <= 1
        O_metric = 1.0;
    else
        mask_upper = triu(true(Q), 1);
        r_off_diag = Phi(mask_upper);
        
        num_pairs = (Q * (Q - 1)) / 2;
        rms_off_diag = sqrt(sum(r_off_diag.^2) / num_pairs);
        O_metric = max(0, min(1, 1 - rms_off_diag)); % Bounded [0, 1]
    end
    
    %% 7. Assemble Output Report Struct
    report = struct();
    
    report.CoreMetrics = struct(...
        'Index', I_index, ...
        'ECV', ECV, ...
        'GS_w', GS_w ...
    );
    report.OrthogonalityCheck = O_metric;
    
    Factor = (1:Q)';
    Weight = w_active;
    Sparseness = gs_vec;
    report.FactorDetails = table(Factor, Weight, Sparseness);
    
    w_Structure = w_str;
    w_Orthogonalized = w_orth;
    w_Direct = w_dir;
    report.WeightSchemes = table(Factor, w_Orthogonalized, w_Structure, w_Direct);
    
    report.Config = struct(...
        'p1', p1, ...
        'p2', p2, ...
        'ActiveWeightScheme', weight_scheme, ...
        'Method', method_label, ...
        'Rotation', rot_label, ...
        'NumVariables', J, ...
        'NumFactors', Q, ...
        'WasPSDProjected', was_projected ...
    );
end

%% Helper Functions
function val = get_opt(opts, field, default_val)
    if isstruct(opts) && isfield(opts, field) && ~isempty(opts.(field))
        val = opts.(field);
    else
        val = default_val;
    end
end

function [w_str, w_orth, w_dir] = compute_all_weights(Lambda, Phi, Phi_sqrt)
    % Direct (Pattern SS)
    ss_cols = sum(Lambda.^2, 1)';
    w_dir = ss_cols / max(sum(ss_cols), eps);
    
    % Structure (Lambda * Phi)
    H = Lambda * Phi;
    sq_H_cols = sum(H.^2, 1)';
    w_str = sq_H_cols / max(sum(sq_H_cols), eps);
    
    % Orthogonalized (Lambda * Phi^(1/2))
    Gamma = Lambda * Phi_sqrt;
    sq_Gamma_cols = sum(Gamma.^2, 1)';
    w_orth = sq_Gamma_cols / max(sum(sq_Gamma_cols), eps);
end