function poses = loadCameraPoses(im_stamps, gt_poses, lcalib, cparam)
im_stamps = im_stamps(:);

x_vb = [0; 0; 0; -pi; 0; 0];
x_vs = ssc_head2tail (x_vb, lcalib.x_vs);

[gt_i, im_i] = findcind(gt_poses.utime, im_stamps);
x_wv = [gt_poses.pos(gt_i, :)'; gt_poses.rph(gt_i, :)'];
x_wv_dot = [gt_poses.vel(gt_i, :)'; gt_poses.rotation_rate(gt_i, :)'];

dt = im_stamps(im_i) - gt_poses.utime(gt_i);
dt = dt * 1e-6;
x_wv = x_wv + x_wv_dot * dt;
x_ws = ssc_head2tail(x_wv, x_vs);

poses(1:3, 1:3, :) = rotxyz(x_ws(4,:), x_ws(5,:), x_ws(6,:));     %[3 x 3 x Nshots]
poses(1:3, 4, :) = reshape(x_ws(1:3, :), 3, 1, []);           %[3 x 1 x Nshots]
poses(4, :, :) = 0; poses(4, 4, :) = 1;

cextrs = inv([cparam.R, cparam.t; 0, 0, 0, 1]);
for i = 1:size(poses, 3)
    poses(:, :, i) = poses(:, :, i) * cextrs;
end