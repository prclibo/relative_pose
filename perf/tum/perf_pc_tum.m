addpath('~/workspace/relative_pose/build/matlab');
addpath('~/workspace/relative_pose/build/Matlab');
addpath('../utils');

clear, clc

data_dir = '~/data/tum/rgbd_dataset_freiburg2_pioneer_slam';
data_dir = '~/data/tum/rgbd_dataset_freiburg2_pioneer_slam2';
data_dir = '~/data/tum/rgbd_dataset_freiburg2_pioneer_slam3';
data_dir = '~/data/tum/rgbd_dataset_freiburg2_pioneer_360';
gt_path = fullfile(data_dir, 'groundtruth.txt');
rgb_list_path = fullfile(data_dir, 'rgb.txt');
rgb_dir = fullfile(data_dir, 'rgb');

fid_rgb_list = fopen(rgb_list_path, 'r');
rgb_list = textscan(fid_rgb_list, '%f %s', 'CommentStyle','#');

rgb_poses = interpolatePoses(gt_path, rgb_list{1}, rgb_list{1}(1));

% https://vision.in.tum.de/data/datasets/rgbd-dataset/file_formats
fx = 525.0;  % focal length x
fy = 525.0;  % focal length y
cx = 319.5;  % optical center x
cy = 239.5;  % optical center y
thresh = 1 / fx;

consec_interv = 1;
start = 50;
rgb_list{1} = rgb_list{1}(start:consec_interv:end);
rgb_list{2} = rgb_list{2}(start:consec_interv:end);
rgb_poses = rgb_poses(start:consec_interv:end);

t_err_5p = nan(size(rgb_list{1}));
t_err_4pst0 = nan(size(rgb_list{1}));
rel_5p = cell(size(rgb_list{1}));
rel_4pst0 = cell(size(rgb_list{1}));
time_5p = nan(size(rgb_list{1}));
time_4pst0 = nan(size(rgb_list{1}));

min_move = 0.1;
gt_move = 0;
for i = 1:numel(rgb_list{1})
    rng(3);
%     clc;
    if i > 1
        offseted = rgb_poses{i};
        offseted_prev = rgb_poses{prev_i};
        offseted(1:3, 4) = offseted(1:3, 4) - offseted_prev(1:3, 4);
        offseted_prev(1:3, 4) = 0;
        gt_rel = offseted_prev \ offseted;
        gt_nt = normc(gt_rel(1:3, 4));
        gt_move = norm(gt_rel(1:3, 4));
    end
    if i == 1 || gt_move > min_move
        curr_path = fullfile(data_dir, rgb_list{2}{i});
        curr_im = rgb2gray(imread(curr_path));
        curr_im = histeq(curr_im, 255);
        curr_points = detectSURFFeatures(curr_im, 'metricthreshold', 300);
        [curr_feat, curr_points] = extractFeatures(curr_im, curr_points);
        if i > 1% && gt_move > min_move
            pairs = matchFeatures(curr_feat, prev_feat,...
                'method', 'Exhaustive', 'maxratio', 0.6);
            mpoints1 = curr_points(pairs(:, 1));
            mpoints2 = prev_points(pairs(:, 2));
            rays1 = ([mpoints1(:).Location] - [cx, cy]) ./ [fx, fy];
            rays2 = ([mpoints2(:).Location] - [cx, cy]) ./ [fx, fy];
            rays1(:, 3) = 1; rays2(:, 3) = 1;
            
            tic, [E_5p, mask_5p] = estimateRelativePose_PC5P_LiH(...
                rays1, rays2, 0.999, thresh);
            time_5p(i) = toc();
            if ~isempty(E_5p)
                pose1 = recoverRelativePose(E_5p, 'rays1', rays1(logical(mask_5p), :), 'rays2', rays2(logical(mask_5p), :));
                if (isfield(pose1, 'R'))
                    t_err_5p(i) = acosd(dot(gt_nt, pose1.t));
                    rel_5p{i} = [pose1.R, pose1.t * gt_move; 0, 0, 0, 1];
                    disp([pose1.R, normc(pose1.t)]);
                end
            end
            
            tic, [E_4pst0, mask_4pst0] = estimateRelativePose_PC4PST0_NullE(...
                rays1, rays2, 0.999, thresh);
            time_4pst0(i) = toc();
            if ~isempty(E_4pst0)
                pose2 = recoverRelativePose(E_4pst0, 'rays1', rays1(logical(mask_4pst0), :), 'rays2', rays2(logical(mask_4pst0), :), 'zeroscrewtransl', true);
                if isfield(pose2, 'R')
                    t_err_4pst0(i) = acosd(dot(gt_nt, pose2.t));
                    rel_4pst0{i} = [pose2.R, pose2.t * gt_move; 0, 0, 0, 1];
                    disp([pose2.R, normc(pose2.t)]);
                end
            end
            disp([sum(mask_4pst0), sum(mask_5p)]);
            disp([gt_rel(1:3, 1:3), normc(gt_rel(1:3, 4)), gt_rel(1:3, 4)]);
            disp('---');
        end
        [prev_points, prev_feat, prev_im] = deal(curr_points, curr_feat, curr_im);
        prev_i = i;
    end
    
%     if i > 300; break; end
end
