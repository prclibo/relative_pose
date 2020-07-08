function pc_sample_data_zp(prime, N, dzero, vars)

assert(~dzero);
rng(3);

sqr = mod((0:prime - 1).^2, prime);
[ii, jj] = meshgrid(1:min(1000, prime), 1:min(1000, prime));
sum_sqr = mod(sqr(ii) + sqr(jj), prime);
[ii, jj] = find(sum_sqr == 1);

assert(numel(ii) > 2);
selected = randi(numel(ii), [1, 3]);

unit = 1;
for i = selected
    unit = [unit(1:end - 1), unit(end) .* [ii(i) - 1, jj(i) - 1]];
    unit = mod(unit, prime);
end
assert(mod(sum(unit.^2), prime) == 1);

unit = sym(unit);
u1 = unit(1);
u2 = unit(2);
u3 = unit(3);
u = [u1; u2; u3];
s = unit(4);
%   R = eval(subs(R, symvar(R), symvar(R)));
R = 2 * (u * transpose(u) - s * skew(u)) + (s * s - transpose(u) * u) * eye(3);
R = mod(R, prime);

t = randi(prime, [3, 1]) - 1; t = sym(t);
Q = randi(prime, [3, N]) - 1; Q = sym(Q);
E = skew(t) * R; E = mod(E, prime);

QQ = R * Q + repmat(t, [1, N]); QQ = mod(QQ, prime);
q = Q; qq = QQ;

[~, invs, ~] = gcd(double(s), prime);
invs = mod(sym(invs), prime);
v1 = mod(u1 * invs, prime);
v2 = mod(u2 * invs, prime);
v3 = mod(u3 * invs, prime);

for r = 1:3
    for c = 1:N
        cmdq = sprintf('q%d%d = sym(%s); ', c, r, char(q(r, c)));
        eval(cmdq);
        cmdqq = sprintf('qq%d%d = sym(%s); ', c, r, char(qq(r, c)));
        eval(cmdqq);
    end
end

names = arrayfun(@(x) char(x), vars, 'UniformOutput', false);
for name = names
    cmd = sprintf('assignin(''base'', ''%s'', %s);', name{1}, name{1});
    eval(cmd);
end
