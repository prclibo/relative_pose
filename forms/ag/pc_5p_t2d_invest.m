% clear, clc

syms u1 u2 u3 s
u = [u1; u2; u3];
unknown_vars = symvar([u; s]);
R = 2 * (u * transpose(u) - s * skew(u)) + (s * s - transpose(u) * u) * eye(3);

clear, clc
v = sym('v%d', [3, 1]);
unknown_vars = symvar(v);

R = 2 * (v * transpose(v) - skew(v)) + (1 - transpose(v) * v) * eye(3);

q = transpose(sym('q%d%d', [5, 3]));
qq = transpose(sym('qq%d%d', [5, 3]));

indices = combinator(5,3,'c');
eqs1 = sym([]);
eqs2 = sym([]);
temp = {};

for index = 1:size(indices, 1)
    i = indices(index, 1);
    j = indices(index, 2);
    k = indices(index, 3);
    F = simpleFijk(q, qq, R, i, j, k);
    
%     G = [R * q(:, i), -qq(:, i), R * q(:, j), -qq(:, j), sym(zeros(3, 2));
%          R * q(:, i), -qq(:, i), sym(zeros(3, 2)), R * q(:, k), -qq(:, k)];
    G = [- transpose(qq(:, i)) * skew(R * q(:, i));
         - transpose(qq(:, j)) * skew(R * q(:, j));
         - transpose(qq(:, k)) * skew(R * q(:, k))];
    
    eqs2(index) = det(G);
    temp{index} = F;
    eqs1(index) = det(F);
    
end

for index = 1:10
    fac = factor(eqs2(index), unknown_vars);
    [c1, t1] = coeffs(eqs1(index), unknown_vars);
    [c2, t2] = coeffs(fac(2), unknown_vars);
    disp(logical(c1 + c2 == 0 & c1 == c2));
    eqs2(index) = fac(2);
end
% 
% eqs3 = transpose(u) * u + s^2 - 1;

%% Verification
% eqs = [eqs1(:); eqs2(:); eqs3(:)];
eqs = eqs2; % subs(eqs1, s, sym(1));

all_vars = [symvar(eqs), str2sym_('[Q, QQ]')];
pc_sample_data(5, false, all_vars);

eqs_ = inf(size(eqs));
for i = 1:numel(eqs)
    fprintf('%d, ', i);
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

