function rotations_struct = manyrot(A, methods, num_starts)
% MANYROT Orchestrator for running multiple factor rotation methods 
%          with timed console progress logging and multi-start details.
%
% Syntax:
%   rotations_struct = manyrot(A)
%   rotations_struct = manyrot(A, 'all')
%   rotations_struct = manyrot(A, {'varimax', 'oblimin', 'promax'}, 50)
    if nargin < 2 || isempty(methods)
        methods = 'all';
    end
    if nargin < 3 || isempty(num_starts)
        num_starts = 25; % Default multi-start iterations
    end
    % Standard benchmark methods in categorized order:
    % 1. Orthogonal: equamax, parsimax, quartimax, varimax
    % 2. Oblique:    oblimin, bentler, promax
    % 3. Modern/Entropy: geomin, infomax, mccammon
    all_methods = { ...
        'equamax', 'parsimax', 'quartimax', 'varimax', ... 
        'promax','oblimin', 'bentler', ...                
        'geomin', 'infomax', 'mccammon' ...  
    };
    
    % Methods that utilize multi-start random initialization in GPA
    multi_start_methods = {'oblimin', 'bentler', 'geomin', 'infomax', 'mccammon'};
    
    if ischar(methods) || isstring(methods)
        if strcmpi(methods, 'all')
            methods_to_run = all_methods;
        else
            methods_to_run = {char(methods)};
        end
    elseif iscell(methods)
        methods_to_run = methods;
    else
        error('Input "methods" must be ''all'', a character vector, or a cell array of strings.');
    end
    
    total_methods = length(methods_to_run);
    rotations_struct = struct();
    
    % Print Header
    fprintf('\n=================================================================\n');
    fprintf('  EXECUTING ROTATIONS (%d methods requested)\n', total_methods);
    fprintf('=================================================================\n');
    
    total_timer = tic; % Start overall benchmark timer
    
    for i = 1:total_methods
        m_name = strtrim(lower(char(methods_to_run{i})));
        
        % Generate a safe struct field name
        field_name = matlab.lang.makeValidName(m_name);
        
        % Determine if this method uses multi-start and format label
        if ismember(m_name, multi_start_methods)
            method_label = sprintf('%s (%d starts)', m_name, num_starts);
        else
            method_label = m_name;
        end
        
        % Print progress status (padded for neat alignment)
        fprintf('  [%2d/%2d] Running %-24s ... ', i, total_methods, method_label);
        
        method_timer = tic; % Start timer for this specific method
        
        try
            % Execute rotation
            [Lambda_rot, T_matlab] = rotate(A, m_name, num_starts);
            
            elapsed_time = toc(method_timer);
            
            % Save outputs and metadata
            rotations_struct.(field_name).Lambda     = Lambda_rot;
            rotations_struct.(field_name).T          = T_matlab;
            rotations_struct.(field_name).runtime_s  = elapsed_time;
            rotations_struct.(field_name).status     = 'SUCCESS';
            
            % Print completion time
            fprintf('Done (%6.3fs)\n', elapsed_time);
            
        catch ME
            elapsed_time = toc(method_timer);
            
            rotations_struct.(field_name).Lambda     = [];
            rotations_struct.(field_name).T          = [];
            rotations_struct.(field_name).runtime_s  = elapsed_time;
            rotations_struct.(field_name).status     = 'FAILED';
            rotations_struct.(field_name).error      = ME.message;
            
            % Handle error cleanly on the same line
            fprintf('FAILED (%6.3fs)\n', elapsed_time);
            warning('Method %s failed: %s', m_name, ME.message);
        end
    end
    
    total_elapsed = toc(total_timer);
    fprintf('=================================================================\n');
    fprintf('All requested rotations completed in %6.2f seconds.\n\n', total_elapsed);
end
