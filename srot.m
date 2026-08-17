function [Lambda, T, f_opt] = srot(A, method_name, num_starts)
% SROT Oblique factor rotation bridge to R's GPArotation package.
%
% Syntax:
%   [Lambda, T] = srot(A, method_name)
%   [Lambda, T] = srot(A, method_name, num_starts)
%   [Lambda, T, f_opt] = srot(...)
%
% Inputs:
%   A           - Unrotated factor loading matrix (Variables x Factors)
%   method_name - Rotation criterion (e.g., 'geomin', 'infomax', 'bentler', 'mccammon')
%   num_starts  - Number of random starts for GPA optimization (default: 25)
%
% Outputs:
%   Lambda      - Primary factor loadings (Variables x Factors)
%   T           - Primary factor transformation matrix (Factors x Factors)
%   f_opt       - Final criterion value at minimum

    %% 1. Input Validation
    narginchk(2, 3);
    validateattributes(A, {'numeric'}, {'2d', 'real', 'finite', 'nonempty'}, mfilename, 'A', 1);
    
    [~, k] = size(A);
    if k < 2
        warning('srot:SingleFactor', 'Matrix A has fewer than 2 factors. Rotation skipped.');
        Lambda = A; T = eye(k); f_opt = 0;
        return;
    end
    
    if nargin < 3 || isempty(num_starts)
        num_starts = 25;
    end
    
    method_str = char(upper(string(method_name)));

    %% 2. Kaiser Normalization
    h2 = sum(A.^2, 2);
    h = sqrt(h2);
    h(h < 1e-10) = 1e-10; % Avoid division by zero
    A_norm = A ./ h;

    %% 3. Setup File Exchange Paths
    user_docs = fileparts(userpath);
    r_folder  = fullfile(user_docs, 'rfiles', 'matlab_bridge');
    
    if ~exist(r_folder, 'dir')
        mkdir(r_folder);
    end
    
    mat_in   = fullfile(r_folder, 'temp_A_norm.mat');
    mat_out  = fullfile(r_folder, 'temp_rot_out.mat');
    r_script = fullfile(r_folder, 'srot.R'); % Corrected script name
    
    if ~exist(r_script, 'file')
        error('srot:MissingRScript', ...
            'R rotation script not found at:\n%s\nPlease ensure srot.R is saved there.', r_script);
    end
    
    if exist(mat_out, 'file'), delete(mat_out); end
    save(mat_in, 'A_norm', 'method_str', 'num_starts', '-v7');

    %% 4. Locate Rscript Executable Dynamically
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

    %% 5. Execute R Script
    cmd = sprintf('%s "%s"', rscript_cmd, r_script);
    curr_dir = cd(r_folder);
    [status, cmdout] = system(cmd);
    cd(curr_dir);

    %% 6. Process Returned Matrices (Mathematically Exact R-to-MATLAB Bridge)
    if status == 0 && exist(mat_out, 'file')
        loaded = load(mat_out);
        
        % loaded.T is R's Reference Vector Transformation Matrix (T_R / fit$Th)
        T_ref = loaded.T; 
        f_opt = loaded.f_opt;
        
        % Convert R's Reference Matrix to MATLAB's Primary Transformation Matrix:
        % T_primary = (T_ref')^-1
        T = (T_ref') \ eye(k);
        
        % Calculate Primary Factor Loadings (Lambda = A_norm * T * h)
        Lambda_norm = A_norm * T;
        Lambda = Lambda_norm .* h;
        
        % Cleanup temporary exchange files
        if exist(mat_in, 'file'), delete(mat_in); end
        if exist(mat_out, 'file'), delete(mat_out); end
    else
        % Print detailed R log if the command failed
        if ~isempty(cmdout)
            fprintf('[srot Diagnostic Console Output]:\n%s\n', cmdout);
        end
        
        warning('srot:FallbackTriggered', ...
            'R rotation (%s) failed to complete. Returning unrotated loadings with identity matrix.', method_str);
        Lambda = A;
        T = eye(k);
        f_opt = NaN;
    end
end