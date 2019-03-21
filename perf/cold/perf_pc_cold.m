clear, rng(3);

addpath(fullfile(fileparts(mfilename('fullpath')), '/../../build/matlab/'));
addpath(fullfile(fileparts(mfilename('fullpath')), '/../utils'));

data_name = 'seq1_cloudy1';
data_dir = fullfile('~/data/cold', data_name);

im_dir = fullfile(data_dir, 'std_cam');
im_list_path = fullfile(data_dir, 'localization/places.lst');

fid_im_list = fopen(im_list_path, 'r');
im_list = textscan(fid_im_list, '%s %s');

im_width = 640;
im_height = 480;
fov_x = deg2rad(68.9);
fov_y = deg2rad(54.4);
fx = im_width / 2 / tan(fov_x / 2);
fy = im_width / 2 / tan(fov_y / 2);
cx = im_width / 2;
cy = im_height / 2;
kc = [-0.2, 0];
cam_param = cameraParameters(...
    'IntrinsicMatrix', [fx, 0, 0; 0, fy, 0; cx, cy, 1], ...
    'RadialDistortion', kc);
thresh = 1 / fx;

cam_extrs = [0, 0, 1, 0; -1, 0, 0, 0; 0, -1, 0, 0; 0, 0, 0, 1];
cam_extrs = [0, 0, 1, 0; 1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 0, 1];

t_err_5p = nan(size(im_list{1}));
t_err_4pst0 = nan(size(im_list{1}));
t_err_4pra = nan(size(im_list{1}));
t_err_3prast0 = nan(size(im_list{1}));
rel_5p = cell(size(im_list{1}));
rel_4pst0 = cell(size(im_list{1}));
rel_4pra = cell(size(im_list{1}));
rel_3prast0 = cell(size(im_list{1}));
time_5p = nan(size(im_list{1}));
time_4pst0 = nan(size(im_list{1}));
time_4pra = nan(size(im_list{1}));
time_3prast0 = nan(size(im_list{1}));

min_move = 0.1;
gt_move = 0;

for i = 1:numel(im_list{1})
    fprintf('%d / %d\n', i, numel(im_list{1}))
    rng(3);
    gt_data = sscanf(im_list{1}{i}, 't%f_x%f_y%f_a%f.jpeg');
    curr_pose = eye(4);
    curr_pose(1:3, 1:3) = vrrotvec2mat([0, 0, 1, gt_data(4)]);
    curr_pose(1:3, 4) = [gt_data(2:3); 0];
    
    if i > 1
        gt_rel = relativePose(curr_pose, prev_pose);
        gt_rel = cam_extrs \ gt_rel * cam_extrs;
        gt_nt = normc(gt_rel(1:3, 4));
        gt_move = norm(gt_rel(1:3, 4));
        gt_E = skew(gt_nt) * gt_rel(1:3, 1:3);
        gt_E = gt_E / norm(gt_E(:));
    end
    if i == 1 || gt_move > min_move
        curr_path = fullfile(im_dir, im_list{1}{i});
        curr_im = rgb2gray(imread(curr_path));
        curr_im = histeq(curr_im, 255);
        curr_im = undistortImage(curr_im, cam_param);
        curr_points = detectSURFFeatures(curr_im, 'metricthreshold', 300);
        [curr_feat, curr_points] = extractFeatures(curr_im, curr_points);
        if i > 1% && gt_move > min_move
            pairs = matchFeatures(curr_feat, prev_feat,...
                'method', 'Exhaustive', 'maxratio', 0.6);
            if size(pairs, 1) > 50
                mpoints1 = curr_points(pairs(:, 1));
                mpoints2 = prev_points(pairs(:, 2));
                %             figure; showMatchedFeatures(curr_im, prev_im, mpoints1, mpoints2);
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
                    end
                end
            end
        end
        [prev_points, prev_feat, prev_im] = deal(curr_points, curr_feat, curr_im);
        prev_i = i;
        prev_pose = curr_pose;
    end
end

% 
% ims = [];
% curr_im = imread('/Users/li/data/cold/seq1_cloudy1/std_cam/t1152886490.040076_x0.205133_y0.000557_a0.002507.jpeg');
% curr_im = rgb2gray(curr_im);
% for k = -0.3:0.05:-0.1
%     cam_param = cameraParameters(...
%         'IntrinsicMatrix', [fx, 0, 0; 0, fy, 0; cx, cy, 1], ...
%         'RadialDistortion', [k, 0]);    
%     und = undistortImage(curr_im, cam_param);
%     ims = [ims, und];
% end
% imshow(ims);