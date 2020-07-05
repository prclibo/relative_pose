clear, clc

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

all_vars = [symvar(eqs), sym('[Q, QQ, R, t]')];
pc_sample_data_zp(cfg.prime, 5, false, all_vars);

F_instance = subs(F, [q(:); qq(:)], eval([q(:); qq(:)]));

eqs_ = sym(inf(size(eqs)));
cfg.eqinstance = sym(inf(size(eqs)));
for r = 1:size(rows, 1)
    fprintf('Instantiating Eq %d\n', r);
    eq_instance = det(F_instance(rows(r, :), :));
    eqs_(r) = mod(subs(eq_instance, unknown_vars, eval(unknown_vars)), cfg.prime);
    cfg.eqinstance(r) = eq_instance;
end
cfg.eqinstance(end) = eqs(end);

for i = 1:numel(cfg.eqinstance)
    eq_instance = cfg.eqinstance(i);
    [coeff, term] = coeffs(eq_instance, unknown_vars);
    coeff = mod(coeff, cfg.prime);
    eq_instance = sum(coeff .* term);
    cfg.eqinstance(i) = eq_instance;

end

% fprintf('Check evaluation of polys on Zp, %f', double(eqs_));

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

all_vars = [symvar(eqs), sym('[Q, QQ]')];
pc_sample_data(5, false, all_vars);

input_str = [sprintf('%s, ', known{1:end - 1}), known{end}];
output_str = [sprintf('%s_, ', unknown{1:end - 1}), sprintf('%s_, ', unknown{end})];

cmd = sprintf('[%s] = solver_%s(%s)', output_str, sname, input_str);
eval(cmd);