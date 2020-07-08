clear, clc

syms u1 u2 u3 s t1 t2
u = [u1; u2; u3];
R = (2 * s * s - 1) * eye(3) + 2 * (u * transpose(u) - s * skew(u));

q = transpose(sym('q%d%d', [4, 3]));
qq = transpose(sym('qq%d%d', [4, 3]));

t = [t1; t2; 1];

E = skew(t) * R;

eqs = sym(zeros(5, 1));
eqs(1) = transpose(qq(:, 1)) * E * q(:, 1);
eqs(2) = transpose(qq(:, 2)) * E * q(:, 2);
eqs(3) = transpose(qq(:, 3)) * E * q(:, 3);
eqs(4) = transpose(qq(:, 4)) * E * q(:, 4);
eqs(5) = transpose(u) * u + s^2 - 1;

unknown_vars = symvar([u, t]);
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

[res, export] = gbs_CreateCode('li4pt_origin', eqs, known, unknown, kngroups);
