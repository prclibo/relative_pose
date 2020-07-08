clear, clc

% syms u1 u2 u3
% u = [u1; u2; u3];
v = sym('v%d', [3, 1]);
unknown_vars = symvar(v);

% R = (2 * s * s - 1) * eye(3) + 2 * (u * transpose(u) - s * skew(u));
R = 2 * (v * transpose(v) - skew(v)) + (1 - transpose(v) * v) * eye(3);

q = transpose(sym('q%d%d', [4, 3]));
qq = transpose(sym('qq%d%d', [4, 3]));

eqs = sym(zeros(4 * 3, 1));
% eqs = sym(zeros(5, 1));

temp = {};

for index = 1:4
    i = index;
    j = mod(index, 4) + 1;
    k = mod(index + 1, 4) + 1;
    F = sym(zeros(3, 2));
    F(1:2, 1:2) = simpleFijk(q, qq, R, i, j, k);
    F(3, :) = [-transpose(v) * q(:, i), transpose(v) * qq(:, i)];
    temp{index} = F;
    eqs(index * 3 - 2) = det(F([1, 2], :));
    eqs(index * 3 - 1) = det(F([2, 3], :));
    eqs(index * 3) = det(F([3, 1], :));
end
% eqs(end) = transpose(u) * u + s^2 - 1;

%% Verification
cfg = gbs_InitConfig();

all_vars = [symvar(eqs), str2sym('[Q, QQ, R, t]')];
pc_sample_data_zp(cfg.prime, 4, false, all_vars);

cfg.eqinstance = sym(inf(size(eqs)));
for r = 1:numel(eqs)
    fprintf('Instantiating Eq %d\n', r);
    eq_instance = subs(eqs(r), [q(:); qq(:)], eval([q(:); qq(:)]));
    cfg.eqinstance(r) = eq_instance;
end

for i = 1:numel(cfg.eqinstance)
    eq_instance = cfg.eqinstance(i);
    [coeff, term] = coeffs(eq_instance, unknown_vars);
    coeff = mod(coeff, cfg.prime);
    eq_instance = sum(coeff .* term);
    cfg.eqinstance(i) = eq_instance;

end


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


