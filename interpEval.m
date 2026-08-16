% =========================================================================
% MASTER DRIVER SCRIPT: EFA & MANY ROTATIONS (INCL. SPARSE)
% =========================================================================
clear; clc; close all;
rng(123)

% -------------------------------------------------------------------------
%% --- DATASET CONFIGURATION ---
dataset_name = 'BostonData'; 
fprintf('Loading dataset: %s...\n', dataset_name);
eval(dataset_name); % Loads raw matrix X
Q = 3;              % Target factor dimension

proc_opts = struct();
proc_opts.is_ordinal = 0; % 0 = Continuous (Pearson); 1 = Ordinal (Polychoric)

%% STEP 1: PRE-PROCESSING & FACTOR EXTRACTION
fprintf('\n--- Step 1: Processing Data & Extracting Factors via processData ---\n');
[Lambda_unrot, R_smooth, Z, data_info] = processData(X, Q, proc_opts);

fprintf('Variables (J): %d | Observations (N): %d | Factors (Q): %d\n', ...
    data_info.NumVariables, data_info.NumObservations, Q);
fprintf('Normality Check Passed : %s\n', string(data_info.IsGaussian));
fprintf('Chosen Extraction Method: %s\n', data_info.ChosenMethod);

%% STEP 2: EXECUTE MULTI-ROTATION EVALUATION VIA MANYROT
fprintf('\n--- Step 2: Executing Factor Rotations via manyrot ---\n');
selected_methods = 'all';

% Initialize structure with baseline Unrotated solution
rotations_struct = struct();
rotations_struct.UNROTATED.Lambda = Lambda_unrot;
rotations_struct.UNROTATED.T      = eye(Q);

% Execute benchmark rotation suite
rot_exec = manyrot(Lambda_unrot, selected_methods);

% Append rotated matrices to master container
fields = fieldnames(rot_exec);
for f = 1:length(fields)
    rotations_struct.(fields{f}) = rot_exec.(fields{f});
end

%% STEP 3: EVALUATE INTERPRETABILITY PROFILE
fprintf('\n--- Step 3: Evaluating Multi-Metric Interpretability ---\n');
opts_base = struct();
opts_base.p1 = 2;
opts_base.p2 = 4;
opts_base.weights = 'orthogonalized';
opts_base.Method = data_info.ChosenMethod;

method_names = fieldnames(rotations_struct);
num_methods  = length(method_names);

methods_col       = cell(num_methods, 1);
index_col         = zeros(num_methods, 1);
ecv_col           = zeros(num_methods, 1);
gsw_col           = zeros(num_methods, 1);
orthogonality_col = zeros(num_methods, 1);
all_reports       = struct();

for i = 1:num_methods
    method_key = method_names{i};
    T_method   = rotations_struct.(method_key).T;
    
    Lambda_rot = Lambda_unrot * T_method;
    Phi_method = inv(T_method' * T_method);
    
    opts_eval = opts_base;
    opts_eval.Rotation = method_key;
    
    report = interp(Lambda_rot, Phi_method, opts_eval);
    all_reports.(method_key) = report;
    
    methods_col{i}       = method_key;
    index_col(i)         = report.CoreMetrics.Index;
    ecv_col(i)           = report.CoreMetrics.ECV;
    gsw_col(i)           = report.CoreMetrics.GS_w;
    orthogonality_col(i) = report.OrthogonalityCheck;
end

% Display final formatted results table
summary_table = table(methods_col, index_col, ecv_col, gsw_col, orthogonality_col, ...
    'VariableNames', {'Method', 'Index', 'ECV', 'GS_w', 'Orthogonality'});

fprintf('\n=======================================================================\n');
fprintf('        FINAL ROTATION COMPARISON TABLE (%s)\n', dataset_name);
fprintf('=======================================================================\n');
disp(summary_table);