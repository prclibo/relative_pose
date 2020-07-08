% clear, clc

syms u1 u2 u3 s
u = [u1; u2; u3];
unknown_vars = symvar([u; s]);
% R = (2 * s * s - 1) * eye(3) + 2 * (u * transpose(u) - s * skew(u));
R = 2 * (u * transpose(u) - s * skew(u)) + (s * s - transpose(u) * u) * eye(3);

q = transpose(sym('q%d%d', [5, 3]));
qq = transpose(sym('qq%d%d', [5, 3]));

F = [- transpose(qq(:, 1)) * skew(R * q(:, 1));
     - transpose(qq(:, 2)) * skew(R * q(:, 2));
     - transpose(qq(:, 3)) * skew(R * q(:, 3));
     - transpose(qq(:, 4)) * skew(R * q(:, 4));
     - transpose(qq(:, 5)) * skew(R * q(:, 5))];

rows = combinator(5,3,'c');

eqs = sym(zeros(size(rows, 1) + 1, 1));
for r = 1:size(rows, 1)
    eqs(r) = det(F(rows(r, :), :));
end
eqs(end) = transpose(u) * u + s^2 - 1;

%% Verification
cfg = gbs_InitConfig();
% cfg.prime = 100003;

all_vars = [symvar(eqs), str2sym_('[Q, QQ, R, t]')];
pc_sample_data(5, true, all_vars);

% eqs_ = inf(size(eqs));
% for i = 1:size(eqs, 1)
%     vars = symvar(eqs(i));
%     eqs_(i) = eval(subs(eqs(i), vars, vars));
% end
% fprintf('eqs evaluated to be less than %e\n', max(abs(eqs_)));


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
[res, export] = gbs_CreateCode(sname, eqs, known, unknown, kngroups, cfg);

%% Post verification
setpaths;
input_str = [sprintf('%s, ', known{1:end - 1}), known{end}];
output_str = [sprintf('%s_, ', unknown{1:end - 1}), sprintf('%s_, ', unknown{end})];

cmd = sprintf('[%s] = solver_%s(%s)', output_str, sname, input_str);
eval(cmd);