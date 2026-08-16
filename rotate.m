function [Lambda, T] = rotate(A, method_name, num_starts)
% ROTATE Executes a single factor rotation returning pattern loadings and T matrix (MATLAB-style).
%
% Syntax:
%   [Lambda, T] = rotate(A, method_name)
%   [Lambda, T] = rotate(A, method_name, num_starts)
%
% Inputs:
%   A           - (p x k) Unrotated factor loadings matrix.
%   method_name - (Char/String) Rotation method name.
%   num_starts  - (Optional) Number of multi-start random initializations (Default = 25).
%
% Outputs:
%   Lambda - (p x k) Rotated factor pattern matrix.
%   T      - (k x k) Transformation matrix (Lambda = A * T).
    if nargin < 3 || isempty(num_starts)
        num_starts = 25; % Default multi-start random initializations
    end
    % Normalize and sanitize method string
    method_str = strtrim(lower(char(method_name)));
    method_upper = upper(method_str);
    switch method_upper
        % --- Native rotatefactors methods (Orthogonal & Oblique Promax) ---
        case {'VARIMAX', 'QUARTIMAX', 'EQUAMAX', 'PARSIMAX', 'PROMAX'}
            % Format for native rotatefactors (expects lowercase/standard method names)
            native_method = lower(method_str);
            [Lambda, T] = rotatefactors(A, 'Method', native_method);
        % --- Direct Oblimin (GPA) ---
        case {'OBLIMIN', 'OBLIMIN_QUARTIMIN'}
            % Pass gamma = 0 (Quartimin) and num_starts to oblimin
            [Lambda, T] = oblimin(A, 0, num_starts);
        case 'OBLIMIN_BIQUARTIMIN'
            % Pass gamma = 0.5 (Biquartimin)
            [Lambda, T] = oblimin(A, 0.5, num_starts);
        % --- Oblique & Modern GPA Rotations ---
        case {'BENTLER', 'GEOMIN', 'INFOMAX', 'MCCAMMON'}
            % Pass method name and num_starts to srot
            [Lambda, T] = srot(A, method_upper, num_starts);
        otherwise
            error('Unrecognized or unsupported factor rotation method "%s".', method_name);
    end
end