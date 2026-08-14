%% --- Configuration ---
WineData; 
Q_max = 13; % Set maximum factors to evaluate
interp_opts = struct('weights', 'orthogonalized');

%% 1. Correlation Matrix & Dimension Setup
[N, J] = size(X);
if N == J && isequal(X, X') && all(abs(diag(X) - 1) < 1e-5)
    R = X;
    N_obs = 500; % Fallback sample size for PA if correlation matrix is passed
else
    R = corr(X, 'Rows', 'pairwise');
    N_obs = N;
end
Q_max = min(Q_max, J);

%% 2. Eigendecomposition & Parallel Analysis
[V, D] = eig(R);
[d_sorted, idx] = sort(diag(D), 'descend');
d_sorted = max(d_sorted, 0); % Prevent tiny negative eigenvalues due to precision
V_sorted = V(:, idx);
eigenvals_data = d_sorted;

n_sims = 50;
sim_eigs = zeros(J, n_sims);
for s = 1:n_sims
    X_null = randn(N_obs, J);
    R_null = corr(X_null);
    sim_eigs(:, s) = sort(eig(R_null), 'descend');
end
eigenvals_pa = mean(sim_eigs, 2);

%% 3. Factor Extraction & Multi-Metric Interpretability Evaluation
ecv_vec   = zeros(Q_max, 1);
gsw_vec   = zeros(Q_max, 1);
index_vec = zeros(Q_max, 1);

for q = 1:Q_max
    % Extract unrotated pattern loading matrix for q factors
    L_q = V_sorted(:, 1:q) * diag(sqrt(d_sorted(1:q)));
    Phi_q = eye(q); % Default unrotated orthogonality
    
    % Call interp.m cleanly
    try
        rep = interp(L_q, Phi_q, interp_opts);
        ecv_vec(q)   = rep.CoreMetrics.ECV;
        gsw_vec(q)   = rep.CoreMetrics.GS_w;
        index_vec(q) = rep.CoreMetrics.Index;
    catch ME
        warning('interp.m failed for q = %d: %s', q, ME.message);
        ecv_vec(q)   = NaN;
        gsw_vec(q)   = NaN;
        index_vec(q) = NaN;
    end
end

%% 4. Summary Table Output
ProfileTable = table((1:Q_max)', ecv_vec, gsw_vec, index_vec, ...
    'VariableNames', {'q', 'ECV', 'GS_w', 'CompositeIndex'});
disp('--- INTERPRETABILITY PROFILE ACROSS Q ---');
disp(ProfileTable);

%% 5. Scree Plot & Interpretability Profile Figures
figure('Color', [1 1 1], 'Position', [100, 100, 1000, 450]);
x_band = [2 4 4 2];

% Panel 1: Scree Plot & Parallel Analysis
subplot(1, 2, 1);
hold on; grid on;
y_lim1 = [0, max(eigenvals_data) * 1.05];

% Highlight Candidate Range Band (Q = 2..4)
patch(x_band, [y_lim1(1) y_lim1(1) y_lim1(2) y_lim1(2)], [0.2 0.5 0.85], ...
    'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

plot(1:J, eigenvals_data, '-o', 'Color', [0.15 0.35 0.65], ...
    'LineWidth', 1.8, 'MarkerFaceColor', [0.15 0.35 0.65], 'DisplayName', 'Data Eigenvalues');
plot(1:J, eigenvals_pa, '--s', 'Color', [0.85 0.35 0.25], ...
    'LineWidth', 1.5, 'MarkerFaceColor', [0.85 0.35 0.25], 'DisplayName', 'Parallel Analysis (Mean)');
line([0.5, J + 0.5], [1 1], 'Color', [0.5 0.5 0.5], 'LineStyle', ':', ...
    'LineWidth', 1.2, 'DisplayName', 'Kaiser Criterion (\lambda=1)');

title('Scree Plot & Parallel Analysis', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Number of Factors (q)', 'FontSize', 10);
ylabel('Eigenvalues', 'FontSize', 10);
xlim([0.5, min(J, Q_max + 0.5)]);
ylim(y_lim1);
legend('Location', 'northeast');

% Panel 2: Multi-Metric Interpretability Profile across q
subplot(1, 2, 2);
hold on; grid on;

% Highlight Candidate Range Band (Q = 2..4)
patch(x_band, [0 0 1 1], [0.2 0.5 0.85], ...
    'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

plot(1:Q_max, ecv_vec, '-^', 'Color', [0.2 0.6 0.3], 'LineWidth', 1.8, ...
    'MarkerFaceColor', [0.2 0.6 0.3], 'DisplayName', 'ECV');
plot(1:Q_max, gsw_vec, '-s', 'Color', [0.8 0.5 0.1], 'LineWidth', 1.8, ...
    'MarkerFaceColor', [0.8 0.5 0.1], 'DisplayName', 'GS_w');
plot(1:Q_max, index_vec, '-o', 'Color', [0.5 0.2 0.7], 'LineWidth', 2.2, ...
    'MarkerFaceColor', [0.5 0.2 0.7], 'DisplayName', 'Composite Index');

title('Interpretability Profile', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Number of Factors (q)', 'FontSize', 10);
ylabel('Metric Score [0, 1]', 'FontSize', 10);
ylim([0, 1]);
xlim([0.5, Q_max + 0.5]);
legend('Location', 'southeast');