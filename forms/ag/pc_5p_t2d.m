clear, clc

syms u1 u2 u3 s
u = [u1; u2; u3];
unknown_vars = symvar([u; s]);
R = (2 * s * s - 1) * eye(3) + 2 * (u * transpose(u) - s * skew(u));

q = transpose(sym('q%d%d', [5, 3]));
qq = transpose(sym('qq%d%d', [5, 3]));

indices = combinator(5,3,'c');
eqs = sym(zeros(size(indices, 1) + 1, 1));

temp = {};

for index = 1:size(indices, 1)
    i = indices(index, 1);
    j = indices(index, 2);
    k = indices(index, 3);
    F = simpleFijk(q, qq, R, i, j, k);
    
    temp{index} = F;
    eqs(index) = det(F);
end
eqs(end) = transpose(u) * u + s^2 - 1;

%% Verification
all_vars = [symvar(eqs), str2sym_('[Q, QQ]')];
pc_sample_data(5, false, all_vars);

for index = 1:1:size(indices, 1)
    vars = symvar(temp{index});
    val = eval(subs(temp{index}, vars, vars));
    i = indices(index, 1);
    
    lambda = norm(Q(:, i));
    mu = norm(QQ(:, i));
    
    disp([val, val(:, 1) ./ val(:, 2)]);
    disp(val * [lambda; mu]);
end
    

eqs_ = inf(size(eqs));
for i = 1:size(eqs, 1)
    vars = symvar(eqs(i));
    eqs_(i) = eval(subs(eqs(i), vars, vars));
end
fprintf('eqs evaluated to be less than %e\n', max(abs(eqs_)));

%% Solve

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
[res, export] = gbs_CreateCode(sname, eqs, known, unknown, kngroups);

%% Post verification
setpaths;
input_str = [sprintf('%s, ', known{1:end - 1}), known{end}];
output_str = [sprintf('%s_, ', unknown{1:end - 1}), sprintf('%s_, ', unknown{end})];

cmd = sprintf('[%s] = solver_%s(%s)', output_str, sname, input_str);
eval(cmd);

% index = 1; v1 = v1_(index); v2 = v2_(index); v3 = v3_(index);
% vars = symvar(eqs(1));
% eeqq = eval(subs(eqs(1), vars, vars))

