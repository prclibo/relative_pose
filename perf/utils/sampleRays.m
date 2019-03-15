function sample = sampleRays(N, noise_std, outlier_rate, angle_turb, planar_turb, motion)

half_view_range = deg2rad(30);
relative_angle_range = deg2rad(5);
axis_turb = deg2rad(5);
transl_turb = deg2rad(5);

angle = randn() * relative_angle_range;

axis = normc([randn() * axis_turb; 1; randn() * axis_turb]);
% axis = [0; 1; 0];

[u, s] = deal(axis .* sin(angle / 2), cos(angle / 2));
angle = angle + randn() * angle_turb;

R = 2 * (u * transpose(u) - s * skew(u)) + (s * s - transpose(u) * u) * eye(3);

if true % parser.Results.ZeroScrewTransl
    if strcmpi(motion, 'forward')
        t = normc([randn([2, 1]) * transl_turb; 1]);
    elseif strcmpi(motion, 'sideway')
        t = normc([1; randn([2, 1]) * transl_turb]);
    else
        error(['Unknown motion: ', motion]);
    end
    t = normc(t - dot(t, axis) * axis);
    t = t + randn() * planar_turb * axis;
else
    error('byebye');
    t = rand([3, 1]);
end
t = normc(t);
% t = [1; 0; 0];


min_depth = 10;
max_depth = 20;
depth = min_depth + rand([N, 1]) * (max_depth - min_depth);

Q = sampleSphereSurface(N, half_view_range);
Q = bsxfun(@times, Q, depth);
Q = Q';

E = skew(t) * R;

QQ = R * Q + repmat(t, 1, N);

M = floor(outlier_rate * N);
QQ(:, N - M + 1:N) = sampleSphereSurface(M, half_view_range)';

q = normc(Q); qq = normc(QQ);
q = q + randn(size(q)) * noise_std;
qq = qq + randn(size(qq)) * noise_std;

names = who();
sample = struct();
for name = names'
    sample.(name{1}) = eval(name{1});
end

end