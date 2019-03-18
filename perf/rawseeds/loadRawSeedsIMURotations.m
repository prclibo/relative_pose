function [imu_rotations, im_indices] = loadRawSeedsIMURotations(imu_path, im_stamps)

imu_data = load(imu_path);
imu_stamps = imu_data(:, 1);

raw_rotations = cell(numel(imu_stamps), 1);
raw_rotations{1} = eye(3);
for i = 2:numel(imu_stamps)
    dt = imu_stamps(i) - imu_stamps(i - 1);
    angular_velocity = mean(imu_data([i - 1, i], 6:8));
    angular_velocity = imu_data(i, 6:8);
    aaxis = angular_velocity * dt;
    inc = vrrotvec2mat([normr(aaxis), norm(aaxis)]);
    raw_rotations{i} = inc * raw_rotations{i - 1};
end

[i_imu, i_im, t_imu, t_im] = findcind(imu_stamps, im_stamps);

nbor_mask = (t_im - t_imu) < 0.1;

i_imu = i_imu(nbor_mask);
i_im = i_im(nbor_mask);

im_indices = i_im;

imu_rotations = raw_rotations(i_imu);