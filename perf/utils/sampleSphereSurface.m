function rays = sampleSphereSurface(N, view_angle)

u = rand([N, 1]);
z = 1 - u * (1 - cos(view_angle));

theta = rand([N, 1]) * 2 * pi;
rho = sqrt(1 - z.^2);
x = cos(theta) .* rho;
y = sin(theta) .* rho;

rays = [x, y, z];

end

