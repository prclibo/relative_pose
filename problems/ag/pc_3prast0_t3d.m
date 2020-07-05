clear, clc

syms u1 u2 u3 s
u = [u1; u2; u3];
R = (2 * s * s - 1) * eye(3) + 2 * (u * transpose(u) - s * skew(u));

q = transpose(sym('q%d%d', [3, 3]));
qq = transpose(sym('qq%d%d', [3, 3]));

F = [- transpose(qq(:, 1)) * skew(R * q(:, 1));
     - transpose(qq(:, 2)) * skew(R * q(:, 2));
     - transpose(qq(:, 3)) * skew(R * q(:, 3));
     transpose(u)];

rows = combinator(4,3,'c');

eqs = sym(zeros(size(rows, 1) + 1, 1));
for r = 1:size(rows, 1)
    eqs(r) = det(F(rows(r, :), :));
end
eqs(end) = transpose(u) * u + s^2 - 1;

unknown_vars = symvar([u]);
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

[res, export] = gbs_CreateCode('pl4pt_det', eqs, known, unknown, kngroups);
