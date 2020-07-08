function pc_sample_data(N, dzero, vars)

rng(3);

unit = rand([4, 1]);
unit = unit ./ norm(unit);

u1 = unit(1);
u2 = unit(2);
u3 = unit(3);
u = [u1; u2; u3];
s = unit(4);
tau = 4 * s^2 - 1;
%   R = eval(subs(R, symvar(R), symvar(R)));
R = (2 * s * s - 1) * eye(3) + 2 * (u * transpose(u) - s * skew(u));

t = rand([3, 1]);
if dzero
    t = t - (dot(t, u) / norm(u)) * (u / norm(u));
    t = normc(t);
end
Q = rand([3, N]);
E = skew(t) * R;

q = normc(Q);
QQ = R * Q + repmat(t, [1, N]);
qq = normc(QQ);

AE = zeros(N, 9);
for i = 1:N
    AE(i, :) = reshape(qq(:, i) * q(:, i)', 1, []);
end
[~, ~, V] = svd(AE);
Evec = V(:, N + 1:end);
ws = Evec \ E(:);
ws = ws ./ ws(end);

v1 = u1 / s;
v2 = u2 / s;
v3 = u3 / s;

for i = 1:numel(ws) - 1
    cmd = sprintf('w%d = %.32f; ', i, ws(i));
    eval(cmd);
end

Evec = reshape(Evec, 3, 3, 9 - N);

% q = bsxfun(@rdivide, q, q(3, :));
% qq = bsxfun(@rdivide, qq, qq(3, :));

for r = 1:3
    for c = 1:N
        cmdq = sprintf('q%d%d = %.32f; ', c, r, q(r, c));
        eval(cmdq);
        cmdqq = sprintf('qq%d%d = %.32f; ', c, r, qq(r, c));
        eval(cmdqq);
    end
end

for i = 1:9 - N
    for r = 1:3
        for c = 1:3
            cmd = sprintf('NE%d%d%d = %.32f;', i, r, c, Evec(r, c, i));
            eval(cmd);
        end
    end
end

names = arrayfun(@(x) char(x), vars, 'UniformOutput', false);
for name = names
    cmd = sprintf('assignin(''caller'', ''%s'', %s);', name{1}, name{1});
    eval(cmd);
end
