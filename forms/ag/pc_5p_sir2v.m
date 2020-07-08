clear, clc

v = sym('v%d', [3, 1]);
unknown_vars = symvar(v);

R = 2 * (v * transpose(v) - skew(v)) + (1 - transpose(v) * v) * eye(3);

q = transpose(sym('q%d%d', [5, 3]));
qq = transpose(sym('qq%d%d', [5, 3]));

indices = combinator(5,3,'c');
eqs = sym(zeros(size(indices, 1), 1));

temp = {};

for index = 1:size(indices, 1)
    i = indices(index, 1);
    j = indices(index, 2);
    k = indices(index, 3);
    F = simpleFijk(q, qq, R, i, j, k);
    
    temp{index} = F;
    eqs(index) = det(F);
end

%% Verification
cfg = gbs_InitConfig();
% cfg.InstanceGenerator = @gbs_RandomInstanceZpFixed;


all_vars = [symvar(eqs), str2sym('[Q, QQ, R, t]')];
pc_sample_data_zp(cfg.prime, 5, false, all_vars);

eqs_ = sym(inf(size(eqs)));
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


