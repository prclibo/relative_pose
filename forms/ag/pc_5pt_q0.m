clear, clc

syms u1 u2 u3 s
u = [u1; u2; u3];
unknown_vars = symvar([u; s]);

R = (2 * s * s - 1) * eye(3) + 2 * (u * transpose(u) - s * skew(u));

q = transpose(sym('q%d%d', [5, 3]));
qq = transpose(sym('qq%d%d', [5, 3]));

eqs = sym(zeros(5 * 6 + 1, 1));

temp = {};

c42 = combinator(4, 2, 'c');

for index = 1:5
    i = index;
    jk = setdiff(1:5, i);
    for jndex = 1:size(c42, 1)
        j = jk(c42(jndex, 1));
        k = jk(c42(jndex, 2));
        disp([i, j, k]);
        F = Fijk(q, qq, R, i, j, k);
        eqs(index * 6 + jndex - 6) = det(F);
        temp{index * 6 + jndex - 6} = F;
    end
end
eqs(end) = transpose(u) * u + s^2 - 1;

%% Verification
unit = rand([4, 1]);
unit = unit ./ norm(unit);

u1 = unit(1);
u2 = unit(2);
u3 = unit(3);
u = [u1; u2; u3];
s = unit(4);
R = eval(subs(R, symvar(R), symvar(R)));

t = rand([3, 1]);
t = t - (dot(t, u) / norm(u)) * (u / norm(u));
Q = rand([3, 5]);
E = skew(t) * R;

q = normc(Q);
tt = R * Q + repmat(t, [1, 5]);
qq = normc(tt);

for r = 1:3
    for c = 1:5
        cmdq = sprintf('q%d%d = %f; ', c, r, q(r, c));
        eval(cmdq);
        cmdqq = sprintf('qq%d%d = %f; ', c, r, qq(r, c));
        eval(cmdqq);
    end
end

for i = 1:size(temp, 1)
    vars = symvar(temp{i});
    val = eval(subs(temp{i}, vars, vars));
    
    lambda = norm(Q(:, i));
    mu = norm(tt(:, i));
    
    disp([val, val(:, 1) ./ val(:, 2)]);
    disp(val * [lambda; mu]);
end
    

eqs_ = inf(size(eqs));
for i = 1:size(eqs, 1)
    vars = symvar(eqs(i));
    eqs_(i) = eval(subs(eqs(i), vars, vars));
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

[res, export] = gbs_CreateCode('pc_5p_q0', eqs, known, unknown, kngroups);
