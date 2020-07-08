clear all
addpath(fullfile('~/workspace/relative_pose/forms', 'helpers/combinator_update/combinator'));
addpath(fullfile('~/workspace/relative_pose/forms', 'probs'));
cd('~/workspace/gaps/');
setpath;
cd('~/workspace/Automatic-Generator/');
setpaths;

prob_fn = @prob_pc_relpose_3prast0_sir2;
prob_fn = @prob_pc_relpose_3prast0_nullex;
% prob_fn = @prob_pc_relpose_5p_sir2;
% prob_fn = @prob_pc_relpose_5p_nulle;
prob_fn = @prob_pc_relpose_4pra_nulle;
prob_fn = @prob_pc_relpose_4pst0_sir2u;
% prob_fn = @prob_pc_relpose_4pst0_sir2v;
prob = prob_fn();

cfg = gbs_InitConfig();
cfg.exportCode = {'matlab'};
% cfg.RemoveUnnecessary = 'fixed';
% cfg.InstanceGenerator = @gbs_RandomInstanceZpFixed;

%% Solve
[~, ~, ~, cfg.eqinstance] = prob.rand_eq_zp(cfg.prime);

sname = char(prob_fn);
eqs = prob.eqs_sym;
eqs = subs_var(eqs, prob.abbr_subs);
known = sym2charcell(cat_flatten_fields(prob.in_subs));
unknown = sym2charcell(cat_flatten_fields(prob.out_subs));
kngroups = [];
[res, export] = gbs_CreateCode(sname, eqs, known(:).', unknown(:).', kngroups, cfg);

%% Verify
rng(23);
[in_rl, out_rl] = prob.rand_arg_rl();
kwn_rl = prob.unpack_pars(prob.in_subs, in_rl);
unk_rl = prob.unpack_pars(prob.out_subs, out_rl);
assign_fields(kwn_rl);
assign_fields(unk_rl);

input_str = [sprintf('%s, ', known{1:end - 1}), known{end}];
output_str = [sprintf('%s_, ', unknown{1:end - 1}), sprintf('%s_, ', unknown{end})];
cmd = sprintf('[%s] = solver_%s(%s)', output_str, sname, input_str);
eval(cmd);