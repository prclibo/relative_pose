addpath(fullfile(fileparts(mfilename('fullpath')), '/../build/matlab/'));
% clear, clc

rng(23);

N = 1000;
K = 5;
noise_std_px = linspace(0, 1, K);
rot_errs_pc4pst0 = zeros(N, K);
rot_errs_pc5p = zeros(N, K);
for k = 1:K
    curr_rot_errs_pc4pst0 = nan(N, 1);
    curr_rot_errs_pc5p = nan(N, 1);
    disp(k);
    for i = 1:N %[26, 35, 60, 82] % 1:N
        rng(k * N + i);
        
        noise_std = noise_std_px(k) / 800;
        sample = sampleRays(5, noise_std);
        
        rays1 = sample.q(:, 1:4)';
        rays2 = sample.qq(:, 1:4)';
        
        E = estimateRelativePose_PC4PST0_NullE(rays1, rays2, 0.99, 1/300);
        pose = recoverRelativePose(E, 'NearestR', sample.R);
        curr_rot_errs_pc4pst0(i) = diff_or_nan(pose, sample);
        
        rays1 = sample.q(:, 1:5)';
        rays2 = sample.qq(:, 1:5)';
        
        E = estimateRelativePose_PC5P_LiH(rays1, rays2, 0.99, 1/300);
        pose = recoverRelativePose(E, 'NearestR', sample.R);
        curr_rot_errs_pc5p(i) = diff_or_nan(pose, sample);
    end
    rot_errs_pc4pst0(:, k) = curr_rot_errs_pc4pst0;
    rot_errs_pc5p(:, k) = curr_rot_errs_pc5p;
end

quantile(rot_errs_pc5p,0.25,1)
quantile(rot_errs_pc4pst0,0.25,1)
