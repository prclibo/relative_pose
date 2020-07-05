clear, clc

R = sym('R%d%d', [3, 3]);
unknown_vars = symvar(R);

q = transpose(sym('q%d%d', [5, 3]));
qq = transpose(sym('qq%d%d', [5, 3]));

indices = combinator(5,3,'c');
eqs1 = sym(zeros(size(indices, 1), 1));

temp = {};

for index = 1:size(indices, 1)
    i = indices(index, 1);
    j = indices(index, 2);
    k = indices(index, 3);
    F = Fijk(q, qq, R, i, j, k);
    
    temp{index} = F;
    eqs1(index) = det(F);
end

eqs2 = [reshape(transpose(R) * R - eye(3), 9, 1);
        reshape(R * transpose(R) - eye(3), 9, 1)];
eqs3 = [skew(R(:, 1)) * R(:, 2) - R(:, 3);
        skew(R(:, 2)) * R(:, 3) - R(:, 1);
        skew(R(:, 3)) * R(:, 1) - R(:, 2)];
eqs = [eqs1; eqs2; eqs3];

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
