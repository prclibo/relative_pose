function sample = sampleRays(N, noise_std)

angle = rand() * deg2rad(45);

axis = normc(rand([3, 1]));
[u, s] = deal(axis .* sin(angle / 2), cos(angle / 2));

R = 2 * (u * transpose(u) - s * skew(u)) + (s * s - transpose(u) * u) * eye(3);

if true % parser.Results.ZeroScrewTransl
    t = skew(u) * rand([3, 1]);
end

t = normc(t);
t = t + rand([3, 1]) * 0.04;

min_depth = 30;
max_depth = 50;
depth = min_depth + rand([N, 1]) * (max_depth - min_depth);

Q = sampleSphereSurface(N, deg2rad(30));
Q = bsxfun(@times, Q, depth);
Q = Q';

E = skew(t) * R;

QQ = R * Q + t;
q = normc(Q); qq = normc(QQ);
q = q + rand(size(q)) * noise_std;
qq = qq + rand(size(qq)) * noise_std;

AE = zeros([N, 9]);
for i = 1:N
    AE(i, :) = reshape(qq(:, i) * q(:, i)', 1, []);
end
nulle = null(AE);
ws = nulle \ E(:);
ws = ws ./ ws(end);

nulle = reshape(nulle, 3, 3, 9 - N);

names = who();
sample = struct();
for name = names'
    sample.(name{1}) = eval(name{1});
end

end