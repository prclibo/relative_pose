addpath(fullfile(fileparts(mfilename('fullpath')), '/../build/matlab/'));
clear, clc



N = 100;
K = 5;
noise_std_px = linspace(0, 1, K);
[rot_errs_pc4pst0, transl_errs_pc4pst0] = deal(zeros(N, K), zeros(N, K));
[rot_errs_pc5p, transl_errs_pc5p] = deal(zeros(N, K), zeros(N, K));
for k = 1:K
    [curr_rot_errs_pc4pst0, curr_transl_errs_pc4pst0] = deal(nan(N, 1), nan(N, 1));
    [curr_rot_errs_pc5p, curr_transl_errs_pc5p] = deal(nan(N, 1), nan(N, 1));
    disp(k);
    for i = 1:N
        rng(k * N + i);
        
        noise_std = noise_std_px(k) / 600;
        sample = sampleRays(5, noise_std, 0);
        
        rays1 = sample.q(:, 1:4)';
        rays2 = sample.qq(:, 1:4)';
        
        E = estimateRelativePose_PC4PST0_NullE(rays1, rays2, 0.99, 1/300);
        pose = recoverRelativePose(E, 'NearestPose', sample);
        curr_rot_errs_pc4pst0(i) = pose.rot_diff;
        curr_transl_errs_pc4pst0(i) = pose.transl_diff;
        
        rays1 = sample.q(:, 1:5)';
        rays2 = sample.qq(:, 1:5)';
        
        E = estimateRelativePose_PC5P_LiH(rays1, rays2, 0.99, 1/300);
        pose = recoverRelativePose(E, 'NearestPose', sample);
        curr_rot_errs_pc5p(i) = pose.rot_diff;
        curr_transl_errs_pc5p(i) = pose.transl_diff;
    end
    rot_errs_pc4pst0(:, k) = curr_rot_errs_pc4pst0;
    transl_errs_pc4pst0(:, k) = curr_transl_errs_pc4pst0;
    rot_errs_pc5p(:, k) = curr_rot_errs_pc5p;
    transl_errs_pc5p(:, k) = curr_transl_errs_pc5p;
end

quantile(rot_errs_pc5p,0.25,1)
quantile(rot_errs_pc4pst0,0.25,1)
quantile(transl_errs_pc5p,0.25,1)
quantile(transl_errs_pc4pst0,0.25,1)