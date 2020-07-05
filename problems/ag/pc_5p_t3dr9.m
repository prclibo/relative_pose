clear, clc

R = sym('R%d%d', [3, 3]);
unknown_vars = symvar(R);

q = transpose(sym('q%d%d', [5, 3]));
qq = transpose(sym('qq%d%d', [5, 3]));

F = [- transpose(qq(:, 1)) * skew(R * q(:, 1));
     - transpose(qq(:, 2)) * skew(R * q(:, 2));
     - transpose(qq(:, 3)) * skew(R * q(:, 3));
     - transpose(qq(:, 4)) * skew(R * q(:, 4));
     - transpose(qq(:, 5)) * skew(R * q(:, 5))];

rows = combinator(5,3,'c');

indices = combinator(5,3,'c');
eqs1 = sym(zeros(size(indices, 1), 1));

temp = {};
for r = 1:size(rows, 1)
    disp(rows(r, :));
    eqs1(r) = det(F(rows(r, :), :));
end

eqs2 = [reshape(transpose(R) * R - eye(3), 9, 1);
        reshape(R * transpose(R) - eye(3), 9, 1)];
eqs3 = [skew(R(:, 1)) * R(:, 2) - R(:, 3);
        skew(R(:, 2)) * R(:, 3) - R(:, 1);
        skew(R(:, 3)) * R(:, 1) - R(:, 2)];
eqs = [eqs1; eqs2; eqs3];

%% Verification
% all_vars = [symvar([q, qq]), sym('[u1, u2, u3, s]')];
% pc_sample_data(5, true, all_vars);


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

