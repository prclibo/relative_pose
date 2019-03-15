function poses = interpolatePoses(gt_path, pose_stamps, origin_stamp)

fid = fopen(gt_path, 'r');
gt_list = textscan(fid, '%f %f %f %f %f %f %f %f', 'CommentStyle','#');
gt_stamps = gt_list{1};

lower_index = max(find(gt_stamps(:,1)<=min(pose_stamps), 1, ...
    'last')-1,1);
upper_index = find(gt_stamps(:,1)>max(pose_stamps), 1, 'first')+1;

lower_index = max(lower_index, 1);
upper_index = min(upper_index, numel(gt_stamps));

gt_poses = cell(1, upper_index - lower_index + 1);
gt_quats = cell(1, upper_index - lower_index + 1);

pose_stamps = pose_stamps(:);
pose_stamps = [origin_stamp; pose_stamps];
gt_stamps = gt_stamps(:);

for i = 1:upper_index
    gt_quats{i} = [-gt_list{8}(i), gt_list{5}(i), gt_list{6}(i), gt_list{7}(i)];
    gt_poses{i} = eye(4);
    gt_poses{i}(1:3, 1:3) = quat2dcm(gt_quats{i});
    gt_poses{i}(1:3, 4) = [gt_list{2}(i), gt_list{3}(i), gt_list{4}(i)];
end

[pose_lindex, gt_lindex] = find(and(...
    bsxfun(@ge, pose_stamps, gt_stamps'), ...
    circshift(bsxfun(@lt, pose_stamps, gt_stamps'), [0, -1])));

fractions = (pose_stamps(pose_lindex) - gt_stamps(gt_lindex))./...
    (gt_stamps(gt_lindex+1)-gt_stamps(gt_lindex));

poses = {};
poses(pose_lindex) = gt_poses(gt_lindex);

for i = numel(poses):-1:1
    poses{i} = poses{1} \ poses{i};
end
poses = poses(2:end);

end