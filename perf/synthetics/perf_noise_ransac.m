addpath(fullfile(fileparts(mfilename('fullpath')), '/../../build/matlab/'));
addpath(fullfile(fileparts(mfilename('fullpath')), '/../utils/'));
addpath(fullfile(fileparts(mfilename('fullpath')), '/../parfor_progress/'));

clear

N = 1000; K1 = 8; K2 = 8; K3 = 8;

focal = 500;
nstd_ray = linspace(0, 1, K1) / focal;
nstd_angle = linspace(0, deg2rad(1), K2);
nstd_transl = linspace(0, 0.05, K3);

rot_errs_pc3prast0 = nan(N, K1, K2, K3); transl_errs_pc3prast0 = nan(N, K1, K2, K3);
rot_errs_pc4pst0 = nan(N, K1, K2, K3); transl_errs_pc4pst0 = nan(N, K1, K2, K3);
rot_errs_pc4pra = nan(N, K1, K2, K3); transl_errs_pc4pra = nan(N, K1, K2, K3);
rot_errs_pc5p = nan(N, K1, K2, K3); transl_errs_pc5p = nan(N, K1, K2, K3);

tic
parfor_progress(N * K1 * K2 * K3);
for i = 1:(N * K1 * K2 * K3)
    rng(i);
    [n, k1, k2, k3] = ind2sub([N, K1, K2, K3], i);
    sample = sampleRays(100, nstd_ray(k1), 0.3, nstd_angle(k2), nstd_transl(k3), 'forward');
    
    rays1 = sample.q';
    rays2 = sample.qq';
    
    E = estimateRelativePose_PC3PRAST0_T2D(sample.angle, rays1, rays2, 0.99, 1 / focal);
    if ~isempty(E)
        pose = recoverRelativePose(E, 'NearestPose', sample);
        rot_errs_pc3prast0(i) = pose.rot_diff;
        transl_errs_pc3prast0(i) = pose.transl_diff;
    end
    
    if k2 == 1
        E = estimateRelativePose_PC4PST0_NullE(rays1, rays2, 0.99, 1 / focal);
        if ~isempty(E)
            pose = recoverRelativePose(E, 'NearestPose', sample);
            rot_errs_pc4pst0(i) = pose.rot_diff;
            transl_errs_pc4pst0(i) = pose.transl_diff;
        end
    end
    if k3 == 1
        E = estimateRelativePose_PC4PRA(sample.angle, rays1, rays2, 0.99, 1 / focal);
        if ~isempty(E)
            pose = recoverRelativePose(E, 'NearestPose', sample);
            rot_errs_pc4pra(i) = pose.rot_diff;
            transl_errs_pc4pra(i) = pose.transl_diff;
        end
    end
    
    if k2 == 1 && k3 == 1
        E = estimateRelativePose_PC5P_LiH(rays1, rays2, 0.99, 1/focal);
        if ~isempty(E)
            pose = recoverRelativePose(E, 'NearestPose', sample);
            rot_errs_pc5p(i) = pose.rot_diff;
            transl_errs_pc5p(i) = pose.transl_diff;
        end
    end
    parfor_progress;
end
parfor_progress(0);
toc

save('/tmp/noise_ransac_motion.mat');