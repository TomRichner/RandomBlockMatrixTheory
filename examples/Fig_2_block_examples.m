% Fig_2_block_examples.m
%
% What block structure buys you: six spectra that the original column-only
% parameterization of Harris et al. (2023) cannot produce.
%
% Every configuration here is one that tests/test_RMTBlocks_predictions.m already
% validates numerically, so the figure and the test agree by construction.
%
% The payoff panel is (f). With f = [0.9 0.1] and a hyper-variable minority, the block
% Perron root and the naive population-averaged variance disagree by ~2x. Both circles
% are drawn: the eigenvalue cloud fills the solid one (ours) and massively overflows the
% dashed one (the naive average). This is Aljadeff, Stern & Sharpee (2015) Fig. 1(a,b) -
% "the average synaptic gain incorrectly predicts this network to be quiescent" - and it
% is the evidence that block structure does not simply average out.
%
% See also: RMTBlocks, Fig_1_RMT_examples

% clear; close all; clc;  % Commented out for master script compatibility

setup_paths();

if exist('master_save_figs', 'var')
    if strcmp(master_save_figs, 'save_all_figs')
        save_figs = true;
    elseif strcmp(master_save_figs, 'save_no_figs')
        save_figs = false;
    end
end
if ~exist('save_figs', 'var')
    save_figs = false;
end

rng(11) % reproducibility

N = 600;        % number of neurons
alpha = 1/3;    % sparse, matching the regime the SRNN work actually runs in
sq = 1/sqrt(N); % all block statistics are given in units of 1/sqrt(N)

%% Configurations
% mu and sigma are (postsynaptic <- presynaptic) blocks:
%   mu(1,2) = mu_EI = "E receives from I"
% Dale's law is a COLUMN constraint - signs constant down each column of mu.

C = struct('name', {}, 'tex', {}, 'note', {}, 'f', {}, 'mu', {}, 'sig', {});

C(end+1) = struct( ...
    'name', 'Column-uniform (Harris)', ...
    'tex',  'Column-uniform (Harris)', ...
    'note', 'reduces to Eq. (17)/(18)', ...
    'f',   [0.5 0.5], ...
    'mu',  [ 3.0 -4.0 ;  3.0 -4.0 ], ...
    'sig', [ 1.0  1.0 ;  1.0  1.0 ]);

C(end+1) = struct( ...
    'name', 'mu_EE elevated', ...
    'tex',  '\mu_{EE} elevated', ...
    'note', 'rank-2 mean -> a second outlier', ...
    'f',   [0.5 0.5], ...
    'mu',  [ 4.5 -4.0 ;  3.0 -4.0 ], ...
    'sig', [ 1.0  1.0 ;  1.0  1.0 ]);

C(end+1) = struct( ...
    'name', 'sigma_II reduced', ...
    'tex',  '\sigma_{II} reduced', ...
    'note', 'block variance moves the bulk edge', ...
    'f',   [0.5 0.5], ...
    'mu',  [ 3.0 -4.0 ;  3.0 -4.0 ], ...
    'sig', [ 1.0  1.0 ;  1.0  0.3 ]);

C(end+1) = struct( ...
    'name', 'Non-separable sigma', ...
    'tex',  'Non-separable \sigma', ...
    'note', 'outside the Ahmadian LJR ensemble', ...
    'f',   [0.5 0.5], ...
    'mu',  [ 0.0  0.0 ;  0.0  0.0 ], ...
    'sig', [ 1.0  0.4 ;  0.4  1.6 ]);

C(end+1) = struct( ...
    'name', 'D = 3 cell types', ...
    'tex',  'D = 3 cell types', ...
    'note', 'general D, three outliers', ...
    'f',   [0.6 0.2 0.2], ...
    'mu',  [ 2.0 -3.0 -1.0 ;  2.0 -3.0 -1.0 ;  1.0 -2.0 -1.0 ], ...
    'sig', [ 1.0  0.8  0.5 ;  0.6  1.2  0.4 ;  0.9  0.3  1.1 ]);

C(end+1) = struct( ...
    'name', 'Hyper-variable minority', ...
    'tex',  'Hyper-variable minority', ...
    'note', 'naive averaged radius fails', ...
    'f',   [0.9 0.1], ...
    'mu',  [ 0.0  0.0 ;  0.0  0.0 ], ...
    'sig', [ 0.5  0.5 ;  0.5  4.0 ]);

n_panels = numel(C);

%% Build
G = cell(1, n_panels);
for k = 1:n_panels
    r = RMTBlocks(N);
    r.alpha = alpha;
    r.set_types(C(k).f, C(k).mu * sq, C(k).sig * sq);
    r.description = C(k).name;
    r.compute_eigenvalues();
    G{k} = r;
end

%% Console summary
fprintf('\n================= Fig 2: block configurations =================\n');
fprintf('N = %d, alpha = %.4f\n\n', N, alpha);
fprintf('%-26s %8s %8s   %s\n', 'configuration', 'R', 'edge', 'lambda_O');
fprintf('%s\n', repmat('-', 1, 78));
for k = 1:n_panels
    r = G{k};
    lam_str = strjoin(arrayfun(@(z) sprintf('%.2f', real(z)), r.lambda_O(:)', ...
        'UniformOutput', false), ', ');
    fprintf('%-26s %8.4f %8.4f   [%s]\n', C(k).name, r.R, r.bulk_edge(99)/sqrt(0.99), lam_str);
end

% Panel (f): the discriminating comparison
r_f = G{end};
ss = r_f.sigma_s_sq;
f_col = r_f.f(:);
R_naive = sqrt(N * (f_col' * ss * f_col));
fprintf('\nPanel (f) discrimination:\n');
fprintf('  measured bulk edge        : %.4f\n', r_f.bulk_edge(99)/sqrt(0.99));
fprintf('  block Perron root (ours)  : %.4f\n', r_f.R);
fprintf('  naive averaged variance   : %.4f   <-- wrong by %.1fx\n', ...
    R_naive, r_f.R / R_naive);
fprintf('===============================================================\n\n');

%% Figure
f2 = figure('Color', 'white');
set(f2, 'Position', [100, 200, 1200, 720]);
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = gobjects(n_panels, 1);
for k = 1:n_panels
    ax(k) = nexttile;
    r = G{k};
    r.plot_spectrum(ax(k));

    hold(ax(k), 'on');

    % Panel (f) only: overlay the naive population-averaged radius, which is what a
    % single-cell-type reading of Harris Eq. (18) would predict.
    if k == n_panels
        theta = linspace(0, 2*pi, 200);
        plot(ax(k), R_naive*cos(theta), R_naive*sin(theta), '--', ...
            'Color', [0.85 0.33 0.10], 'LineWidth', 2);
        text(ax(k), 0, -R_naive, sprintf('  naive avg = %.2f', R_naive), ...
            'Color', [0.85 0.33 0.10], 'FontSize', 9, ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    end
    hold(ax(k), 'off');

    % Frame each panel around its own spectrum, with room for the outliers
    lam = r.eigenvalues();
    lim_x = max(abs(real(lam))) * 1.15;
    lim_y = max(max(abs(imag(lam))), r.R) * 1.15;
    lim_both = max(lim_x, lim_y);
    xlim(ax(k), [-lim_both, lim_both]);
    ylim(ax(k), [-lim_both, lim_both]);

    lam_str = strjoin(arrayfun(@(z) sprintf('%.1f', real(z)), r.lambda_O(:)', ...
        'UniformOutput', false), ', ');
    title(ax(k), {C(k).tex, ...
        sprintf('\\rm\\fontsize{9}%s', C(k).note), ...
        sprintf('\\rm\\fontsize{9}R = %.2f,  \\lambda_O = [%s]', r.R, lam_str)});
    grid(ax(k), 'on');
end

letters = arrayfun(@(c) sprintf('(%s)', c), 'a':'z', 'UniformOutput', false);
AddLetters2Plots(num2cell(ax), letters, 'FontSize', 16, 'FontWeight', 'normal', ...
    'HShift', +0.01, 'VShift', -0.03);

drawnow;

%% Supplementary: the W matrices, where the block structure is directly visible
f3 = figure('Color', 'white');
set(f3, 'Position', [150, 150, 1200, 720]);
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax3 = gobjects(n_panels, 1);
for k = 1:n_panels
    ax3(k) = nexttile;
    G{k}.plot_W(ax3(k));
    title(ax3(k), C(k).tex, 'Color', 'black');
end

drawnow;

%% Save
if save_figs
    fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figs');
    save_some_figs_to_folder_2(fig_dir, 'block_examples', [f2.Number, f3.Number], ...
        {'fig', 'svg', 'png'});
end
