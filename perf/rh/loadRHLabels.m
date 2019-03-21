function [gt_poses, im_ids] = loadRHLabels(label_path, sensor_label)

fid = fopen(label_path, 'r');

label_data = textscan(fid, '%d %s %f %f %f %f %f %f %f', ...
    'CommentStyle', '#');

sensor_indices = find(strcmpi(label_data{2}, sensor_label));
im_ids = label_data{1}(sensor_indices);


[~, si] = sort(label_data{9});
sl = label_data{2}(si);
x = label_data{3}(si);
y = label_data{4}(si);
z = label_data{5}(si);
yaws = label_data{6}(si);
pitchs = label_data{7}(si);
rolls = label_data{8}(si);

% jj = sensor_indices(1);
% yaw = yaws(jj);pitch = pitchs(jj); roll = rolls(jj);
% 
% R0 = [cos(yaw) * cos(pitch), cos(yaw)*sin(pitch)*sin(roll) - sin(yaw)*cos(roll), cos(yaw)*sin(pitch)*cos(roll) + sin(yaw)*sin(roll);
%            sin(yaw) * cos(pitch), sin(yaw)*sin(pitch)*sin(roll) + cos(yaw)*cos(roll), sin(yaw)*sin(pitch)*cos(roll) - cos(yaw)*sin(roll);
%            -sin(pitch), cos(pitch)*sin(roll), cos(pitch)*cos(roll)];
% t0 = [x(jj); y(jj); z(jj)];
% jj = sensor_indices(36);
% yaw = yaws(jj);pitch = pitchs(jj); roll = rolls(jj);
% 
% R1 = [cos(yaw) * cos(pitch), cos(yaw)*sin(pitch)*sin(roll) - sin(yaw)*cos(roll), cos(yaw)*sin(pitch)*cos(roll) + sin(yaw)*sin(roll);
%            sin(yaw) * cos(pitch), sin(yaw)*sin(pitch)*sin(roll) + cos(yaw)*cos(roll), sin(yaw)*sin(pitch)*cos(roll) - cos(yaw)*sin(roll);
%            -sin(pitch), cos(pitch)*sin(roll), cos(pitch)*cos(roll)];
% t1 = [x(jj); y(jj); z(jj)];


gt_poses = cell(numel(sensor_indices), 1);
for i = 1:numel(sensor_indices)
    si = sensor_indices(i);
    gt_poses{i} = eye(4);
    gt_poses{i}(1:3, 1:3) = angle2dcm(...
        label_data{6}(si), label_data{7}(si), label_data{8}(si))';
    
    yaw = label_data{6}(si);
    pitch = label_data{7}(si);
    roll = label_data{8}(si);
    R = [cos(yaw) * cos(pitch), cos(yaw)*sin(pitch)*sin(roll) - sin(yaw)*cos(roll), cos(yaw)*sin(pitch)*cos(roll) + sin(yaw)*sin(roll);
           sin(yaw) * cos(pitch), sin(yaw)*sin(pitch)*sin(roll) + cos(yaw)*cos(roll), sin(yaw)*sin(pitch)*cos(roll) - cos(yaw)*sin(roll);
           -sin(pitch), cos(pitch)*sin(roll), cos(pitch)*cos(roll)];
%     disp([R, gt_poses{i}(1:3, 1:3)]);
       

    
%     gt_poses{i}(1:3, 1:3) = gt_poses{i}(1:3, 1:3) * ...
%         [0, 0, 1; 1, 0, 0; 0, 1, 0];
    gt_poses{i}(1:3, 4) = [...
        label_data{3}(si); label_data{4}(si); label_data{5}(si)];
end

end