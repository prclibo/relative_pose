addpath(fullfile(fileparts(mfilename('fullpath')), '/../../matlab/'));
addpath(fullfile(fileparts(mfilename('fullpath')), '/../../build/matlab/'));

clear

N = 100; K1 = 3; K2 = 2; K3 = 2;

focal = 500;
nstd_ray = linspace(0, 1, K1) / focal;
nstd_angle = linspace(0, deg2rad(5), K2);
nstd_transl = linspace(0, 0.05, K3);

nroots_pc3prast0 = zeros(N, K1, K2, K3);
nroots_pc4pst0 = zeros(N, K1, K2, K3);
nroots_pc4pra = zeros(N, K1, K2, K3);
nroots_pc5p = zeros(N, K1, K2, K3);

for i = 1:(N * K1 * K2 * K3)
    fprintf('%d / %d\n', i, N * K1 * K2 * K3);
    rng(i);
    [n, k1, k2, k3] = ind2sub([N, K1, K2, K3], i);
    sample = sampleRays(5, nstd_ray(k1), 0.3, nstd_angle(k2), nstd_transl(k3), 'sideway');
    
    rays1 = sample.q';
    rays2 = sample.qq';
    
    E = estimateRelativePose_PC3PRAST0_T2D(sample.angle, rays1(1:3, :), rays2(1:3, :), 0.99, 1 / focal / 2);
    nroots_pc3prast0(i) = size(E, 1) / 3;
    
%     if k2 == 1
        E = estimateRelativePose_PC4PST0_NullE(rays1(1:4, :), rays2(1:4, :), 0.99, 1 / focal / 2);
        nroots_pc4pst0(i) = size(E, 1) / 3;
%     end
%     if k3 == 1
        E = estimateRelativePose_PC4PRA(sample.angle, rays1(1:4, :), rays2(1:4, :), 0.99, 1 / focal / 2);
        nroots_pc4pra(i) = size(E, 1) / 3;
%     end
    
%     if k2 == 1 && k3 == 1
        E = estimateRelativePose_PC5P_LiH(rays1, rays2, 0.99, 1/focal / 2);
        nroots_pc5p(i) = size(E, 1) / 3;
%     end
end
toc
