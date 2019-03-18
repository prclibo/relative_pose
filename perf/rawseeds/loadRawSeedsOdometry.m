function [odo_poses, im_indices] = loadRawSeedsOdometry(odo_path, im_stamps)

odo_data = load(odo_path);

odo_stamps = odo_data(:, 1);

[i_odo, i_im, t_odo, t_im] = findcind(odo_stamps, im_stamps);

nbor_mask = (t_im - t_odo) < 0.1;

i_odo = i_odo(nbor_mask);
i_im = i_im(nbor_mask);

im_indices = i_im;

odo_poses = cell(numel(im_indices), 1);
for i = 1:numel(i_odo)
    odo_poses{i} = eye(4);
    odo_poses{i}(1:3, 1:3) = vrrotvec2mat([0, 0, 1, odo_data(i_odo(i), 7)]);
    odo_poses{i}(1:3, 4) = [odo_data(i_odo(i), 5:6)'; 0];
end