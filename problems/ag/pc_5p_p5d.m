clear, clc

syms v1 v2 v3 t1 t2
v = [v1; v2; v3];
unknown_vars = symvar([v; t1; t2]);
R = 2 * (v * transpose(v) - skew(v)) + (1 - transpose(v) * v) * eye(3);

q = transpose(sym('q%d%d', [5, 3]));
qq = transpose(sym('qq%d%d', [5, 3]));

t = [t1; t2; 1];

E = skew(t) * R;

eqs = sym(zeros(5, 1));
eqs(1) = transpose(qq(:, 1)) * E * q(:, 1);
eqs(2) = transpose(qq(:, 2)) * E * q(:, 2);
eqs(3) = transpose(qq(:, 3)) * E * q(:, 3);
eqs(4) = transpose(qq(:, 4)) * E * q(:, 4);
eqs(5) = transpose(qq(:, 5)) * E * q(:, 5);

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
