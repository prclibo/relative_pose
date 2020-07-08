clear all
addpath(fullfile('~/workspace/relative_pose/forms', 'helpers/combinator_update/combinator'));
addpath(fullfile('~/workspace/relative_pose/forms', 'probs'));
cd('~/workspace/gaps/');
setpath;
cd('~/workspace/Automatic-Generator/');
setpaths;

opt = default_options();
opt.M2_path = '/Users/li/workspace/Macaulay2-1.13/bin/M2';
opt.M2_path = 'M2';
opt.optimize_coefficients = true;
opt.remove_extra_columns = false;
opt.find_upper_trianglar = false;
% opt.sparse_template = true;
% opt.cg_language = 'cpp';
% opt.eigen_solver = 'eigs_only';
opt.eigen_solver = 'default';
opt.use_sym = true;
opt.cg_eigen_dir = '/usr/local/include/eigen3';

% opt.find_upper_trianglar = false;

prob_fn = @prob_pc_relpose_4pra_sir2u;
prob_fn = @prob_pc_relpose_3prast0_sir2u;
prob_fn = @prob_pc_relpose_5p_sir2u;
% prob_fn = @prob_pc_relpose_5p_sir3u;
% prob_fn = @prob_pc_relpose_5p_nulle;
% prob_fn = @prob_pc_relpose_3prast0_nullex;
% prob_fn = @prob_pc_relpose_4pra_nulle;
prob_fn = @prob_pc_relpose_4pst0_sir2u;
% prob_fn = @prob_pc_relpose_4pst0_sir2v;
% prob_fn = @prob_pc_relpose_5p_sir2v;
[solv, opt] = generate_solver(prob_fn, opt);
% return;
addpath solvers

solv_fun = str2func(['solver_' solv.name]);
stats = benchmark_solver(solv_fun,solv.prob,500);

figure(1)
clf
hist(log10(stats.all_res),50)
title(sprintf('Mean: %.2f, Median: %.2f\n Mode: %.2f, Time: %.2f ms\n',...
        stats.res_mean,stats.res_median,stats.res_mode,...
        1000*median(stats.time_taken)));
xlabel('log10 residual')
ylabel('freq.')
