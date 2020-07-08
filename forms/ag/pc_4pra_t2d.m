clear, clc

syms u1 u2 u3 s
u = [u1; u2; u3];
R = (2 * s * s - 1) * eye(3) + 2 * (u * transpose(u) - s * skew(u));

q = transpose(sym('q%d%d', [4, 3]));
qq = transpose(sym('qq%d%d', [4, 3]));


% gc_zp_data(4, false, [symvar([q, qq]), sym('s')]);
% q = eval(q);
% qq = eval(qq);
% R = eval(R);

eqs = sym(zeros(5, 1));
eqs(1) = det(simpleFijk(q, qq, R, 2, 3, 4));
eqs(2) = det(simpleFijk(q, qq, R, 3, 4, 1));
eqs(3) = det(simpleFijk(q, qq, R, 4, 1, 2));
eqs(4) = det(simpleFijk(q, qq, R, 1, 2, 3));
eqs(5) = transpose(u) * u + s^2 - 1;

%% Verification
pc_sample_data(4, true, setdiff([symvar([q, qq, R]), transpose(u), s], []));
% u = eval(u);
% q = eval(q);
% qq = eval(qq);
% R = eval(R);

eqs_ = inf(size(eqs));
for i = 1:size(eqs, 1)
    vars = symvar(eqs(i));
    eqs_(i) = eval(subs(eqs(i), vars, vars));
end
fprintf('eqs evaluated to be less than %e\n', max(abs(eqs_)));

%% Solve

unknown_vars = str2sym_('[u1, u2, u3]');
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