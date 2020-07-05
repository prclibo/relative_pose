clear, clc

syms u1 u2 u3 s
u = [u1; u2; u3];
unknown_vars = symvar([u; s]);
% R = (2 * s * s - 1) * eye(3) + 2 * (u * transpose(u) - s * skew(u));
R = 2 * (u * transpose(u) - s * skew(u)) + (s * s - transpose(u) * u) * eye(3);

q = transpose(sym('q%d%d', [4, 3]));
qq = transpose(sym('qq%d%d', [4, 3]));


% gc_zp_data(4, false, [symvar([q, qq]), sym('s')]);
% q = eval(q);
% qq = eval(qq);
% R = eval(R);

F = {};
F{1} = Fijk(q, qq, R, 2, 3, 4);
F{2} = Fijk(q, qq, R, 3, 4, 1);
F{3} = Fijk(q, qq, R, 4, 1, 2);
F{4} = Fijk(q, qq, R, 1, 2, 3);

eqs = sym(zeros(5, 1));
eqs(1) = det(F{1});
eqs(2) = det(F{2});
eqs(3) = det(F{3});
eqs(4) = det(F{4});
eqs(5) = transpose(u) * u + s^2 - 1;

%% Verification
cfg = gbs_InitConfig();

all_vars = [symvar(eqs), sym('[Q, QQ, R, t]')];
pc_sample_data_zp(cfg.prime, 4, false, all_vars);

eqs_ = sym(inf(4, 1));
cfg.eqinstance = sym(inf(size(eqs)));
for r = 1:4
    fprintf('Instantiating Eq %d\n', r);
    F_instance = subs(F{r}, [q(:); qq(:)], eval([q(:); qq(:)]));
    eq_instance = det(F_instance);
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

%% Solve

unknown_vars = sym('[u1, u2, u3]');
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

all_vars = [symvar(eqs), sym('[Q, QQ]')];
pc_sample_data(4, false, all_vars);

input_str = [sprintf('%s, ', known{1:end - 1}), known{end}];
output_str = [sprintf('%s_, ', unknown{1:end - 1}), sprintf('%s_, ', unknown{end})];
cmd = sprintf('[%s] = solver_%s(%s)', output_str, sname, input_str);
eval(cmd);