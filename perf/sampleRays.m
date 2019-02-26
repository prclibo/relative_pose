function sample = sampleRays(N, noise_std, outlier_rate)

half_view_range = deg2rad(30);
relative_angle_range = deg2rad(2);
transl_turb = 0.01;

angle = rand() * relative_angle_range;

axis = normc(rand([3, 1]));
[u, s] = deal(axis .* sin(angle / 2), cos(angle / 2));

R = 2 * (u * transpose(u) - s * skew(u)) + (s * s - transpose(u) * u) * eye(3);

if true % parser.Results.ZeroScrewTransl
    t = skew(u) * rand([3, 1]);
    t = normc(t);
    t = t + rand([3, 1]) * transl_turb;
else
    t = rand([3, 1]);
end
t = normc(t);


min_depth = 30;
max_depth = 50;
depth = min_depth + rand([N, 1]) * (max_depth - min_depth);

Q = sampleSphereSurface(N, half_view_range);
Q = bsxfun(@times, Q, depth);
Q = Q';

E = skew(t) * R;

QQ = R * Q + repmat(t, 1, N);

M = floor(outlier_rate * N);
QQ(:, N - M + 1:N) = sampleSphereSurface(M, half_view_range)';

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