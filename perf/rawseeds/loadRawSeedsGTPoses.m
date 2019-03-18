function [gt_poses, im_indices] = loadRawSeedsGTPoses(gt_path, im_stamps)

gt_data = load(gt_path);
gt_stamps = gt_data(:, 1);

[i_gt, i_im, t_gt, t_im] = findcind(gt_stamps, im_stamps);

nbor_mask = (t_im - t_gt) < 0.1;

i_gt = i_gt(nbor_mask);
i_im = i_im(nbor_mask);

im_indices = i_im;

gt_poses = cell(numel(im_indices), 1);
for i = 1:numel(i_gt)
    gt_poses{i} = eye(4);
    gt_poses{i}(1:3, 1:3) = vrrotvec2mat([0, 0, 1, gt_data(i_gt(i), 4)]);
    gt_poses{i}(1:3, 4) = [gt_data(i_gt(i), 2:3)'; 0];
end
