clear, clc
E1 = sym('NE1%d%d', [3 ,3]);
E2 = sym('NE2%d%d', [3 ,3]);
E3 = sym('NE3%d%d', [3 ,3]);
E4 = sym('NE4%d%d', [3 ,3]);

syms w1 w2 w3
E = w1*E1 + w2*E2 + w3*E3 + E4;

eqs = sym(zeros(10, 1));
eqs(1) = det(E);

Et = transpose(E);
te = 2*(E*Et)*E - trace(E*Et)*E;
eqs(2:10) = te(:);

%% Null verification
all_vars = [symvar(eqs), str2sym_('[Q, QQ]')];
pc_sample_data(5, true, all_vars);

eqs_ = inf(size(eqs));
for i = 1:size(eqs, 1)
    vars = symvar(eqs(i));
    eqs_(i) = eval(eqs(i));
end
fprintf('eqs evaluated to be less than %e\n', max(abs(eqs_)));

%% Solve

unknown_vars = sym('w%d', [1, 3]);

known_vars = setdiff(symvar(eqs), unknown_vars);
known = {};
for var = known_vars
    known = [known {char(var)}];
end
unknown = {};
for var = unknown_vars
    unknown = [unknown {char(var)}];
end

kngroups = [];
sname = mfilename();
cfg = gbs_InitConfig();
cfg.InstanceGenerator = @gbs_RandomInstanceZpFixed;
[res, export] = gbs_CreateCode(sname, eqs, known, unknown, kngroups, cfg);

%% Post verification
setpaths;
input_str = [sprintf('%s, ', known{1:end - 1}), known{end}];
output_str = [sprintf('%s_, ', unknown{1:end - 1}), sprintf('%s_, ', unknown{end})];

cmd = sprintf('[%s] = solver_%s(%s)', output_str, sname, input_str);
eval(cmd);

% fprintf('sampled v: [%f, %f, %f] \n', v1, v2, v3);
% index = 1; v1 = v1_(index); v2 = v2_(index); v3 = v3_(index);
% vars = symvar(eqs(1));
% eeqq = eval(subs(eqs(1), vars, vars))

