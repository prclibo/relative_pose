clear;
rng(3);

% 5.1 & 7.9 \\ 
% 10.3 & 14.3 \\ 
% 15.7 & 22.8 \\ 
% 
% 3.4 & 7.2 \\ 
% 9.3 & 12.0 \\ 
% 16.3 & 16.9 \\ 
% 
% 4.9 & 5.3 \\ 
% 8.5 & 9.9 \\ 
% 16.1 & 23.5 \\ 
% 
% 6.3 & 6.1 \\ 
% 10.7 & 10.7 \\ 
% 17.9 & 17.7 \\ 
% 
% 7.3 & 9.0 \\ 
% 12.2 & 17.4 \\ 
% 21.4 & 28.0 \\ 
% 
% 
% 8.7 & 9.2 \\ 
% 14.3 & 21.5 \\ 
% 22.9 & 39.3 \\ 
% 
% 8.1 & 8.5 \\ 
% 13.3 & 15.4 \\ 
% 22.4 & 24.8 \\ 
% 
% 5.0 & 6.4 \\ 
% 8.6 & 11.7 \\ 
% 13.7 & 21.5 \\ 
% 
% 7.7 & 13.3 \\ 
% 13.4 & 19.1 \\ 
% 30.8 & 26.5 \\ 
% 
% 9.4 & 12.9 \\ 
% 17.1 & 29.8 \\ 
% 34.7 & 46.9 \\ 


% anto-s1
%     6.7267    7.3874
%    11.2065   13.2670
%    20.4855   29.5818
% alma-s1
%     9.3131    8.4759
%    14.6076   14.8671
%    39.7566   40.8541
% pare-s1
%     5.6766    5.9751
%    11.3704   12.3091
%    25.8332   39.6318
% rx2-s1
%     5.6766    5.9751
%    11.3704   12.3091
%    25.8332   39.6318
% sarmis-s1
%     6.4396    7.3851
%    12.6020   16.4660
%    25.2733   53.2055

clear;rng(3);

addpath('~/data/umich_ford/Code/MATLAB/create_ijrr_utils/')
addpath(fullfile(fileparts(mfilename('fullpath')), '/../../build/matlab/'));
addpath(fullfile(fileparts(mfilename('fullpath')), '/../utils'));

choose_data = 1;
choose_cam = 'side';
if choose_data == 1
    data_name = 'anto-s1';
    if strcmp(choose_cam, 'front'); sensor_label = 'RGBD_4';
    elseif strcmp(choose_cam, 'side'); sensor_label = 'RGBD_1'; end
    room_names = {'bathroom1_labelled', 'bedroom1_labelled', 'corridor1_labelled', 'kitchen1_labelled', 'masterroom1_labelled', 'bathroom2_labelled', 'bedroom2_labelled', 'dressingroom1_labelled', 'livingroom1_labelled'};
elseif choose_data == 2
    data_name = 'alma-s1';
    if strcmp(choose_cam, 'front'); sensor_label = 'RGBD_4'; 
    elseif strcmp(choose_cam, 'side'); sensor_label = 'RGBD_1'; end
    room_names = {'bathroom1_labelled', 'bedroom1_labelled', 'kitchen1_labelled', 'livingroom1_labelled', 'masterroom1_labelled' };
elseif choose_data == 3
    data_name = 'pare-s1';
    if strcmp(choose_cam, 'front'); sensor_label = 'RGBD_4'; 
    elseif strcmp(choose_cam, 'side'); sensor_label = 'RGBD_1'; end
    room_names = {'bathroom1_labelled', 'bedroom1_labelled', 'corridor1_labelled', 'kitchen1_labelled', 'livingroom2_labelled', 'bathroom2_labelled', 'bedroom2_labelled', 'hall1_labelled', 'livingroom1_labelled', 'masterroom1_labelled' };
elseif choose_data == 4
    data_name = 'rx2-s1';
    if strcmp(choose_cam, 'front'); sensor_label = 'RGBD_4'; 
    elseif strcmp(choose_cam, 'side'); sensor_label = 'RGBD_1'; end
    room_names = {'bathroom1_labelled', 'bedroom1_labelled', 'kitchen1_labelled', 'livingroom1_labelled' };
elseif choose_data == 5
    data_name = 'sarmis-s1';
    if strcmp(choose_cam, 'front'); sensor_label = 'RGBD_1'; 
    elseif strcmp(choose_cam, 'side'); sensor_label = 'RGBD_2'; end
    room_names = {'bathroom2_1_labelled', 'bedroom1_1_labelled', 'bedroom2_1_labelled', 'bedroom3_1_labelled', 'corridor1_1_labelled', 'kitchen1_1_labelled', 'livingroom1_1_labelled', 'bathroom2_2_labelled', 'bedroom1_2_labelled', 'bedroom2_2_labelled', 'bedroom3_2_labelled', 'corridor1_2_labelled', 'kitchen1_2_labelled', 'livingroom1_2_labelled'};
    
end
data_dir = fullfile('~/data/rh', data_name);

cx = 157.3245865;
cy = 120.0802295;
fx = 572.882768/2;
fy = 542.739980/2;

cam_extrs = [0, 0, 1, 0; -1, 0, 0, 0; 0, -1, 0, 0; 0, 0, 0, 1];
% cam_extrs = [0, -1, 0, 0; 1, 0, 0, 0; 0, 0, 1, 0;0, 0, 0, 1];

thresh = 1 / fx;

t_err_5p = [];
t_err_4pst0 = [];

min_move = 0.1;

for room_i = 1:numel(room_names)
%     if room_i > 5; return; end
    im_dir = fullfile(data_dir, room_names{room_i});
    label_path = [im_dir, '.txt'];
    [gt_poses, im_ids] = loadRHLabels(label_path, sensor_label);
    [fr_poses, fr_ids] = loadRHLabels(label_path, 'RGBD_4');
    curr_t_err_5p = nan(size(im_ids));
    curr_t_err_4pst0 = nan(size(im_ids));
    gt_move = 0;
    prev_i = 0;
    for i = 1:numel(im_ids)
        fprintf('%d / %d\n', i, numel(im_ids))
        rng(3);
        %     clc
        if i > 1
            gt_rel = relativePose(gt_poses{i}, gt_poses{prev_i});
            gt_rel = cam_extrs \ gt_rel * cam_extrs;
            gt_nt = normc(gt_rel(1:3, 4));
            gt_move = norm(gt_rel(1:3, 4));
            gt_E = skew(gt_nt) * gt_rel(1:3, 1:3);
            gt_E = gt_E / norm(gt_E(:));
        end
        if i == 1 || gt_move > min_move
            curr_path = fullfile(im_dir, sprintf('%d_intensity.png', im_ids(i)));
            curr_im = rgb2gray(imread(curr_path));
            curr_im = histeq(curr_im, 255);
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
                    if ~isempty(E_5p)
                        pose1 = recoverRelativePose(E_5p, 'rays1', rays1(logical(mask_5p), :), 'rays2', rays2(logical(mask_5p), :));
                        if (isfield(pose1, 'R'))
                            curr_t_err_5p(i) = acosd(dot(gt_nt, pose1.t));
                        end
                    end
                    tic, [E_4pst0, mask_4pst0] = estimateRelativePose_PC4PST0_NullE(...
                        rays1, rays2, 0.999, thresh);
                    if ~isempty(E_4pst0)
                        pose2 = recoverRelativePose(E_4pst0, 'rays1', rays1(logical(mask_4pst0), :), 'rays2', rays2(logical(mask_4pst0), :), 'zeroscrewtransl', true);
                        if isfield(pose2, 'R')
                            curr_t_err_4pst0(i) = acosd(dot(gt_nt, pose2.t));
                        end
                    end
                    [E_2pot, mask_2pot] = estimateRelativePose_PC2POT(...
                        rays1, rays2, 0.999, thresh);
                    if ~isempty(E_2pot)
                        pose2ot = recoverRelativePose(E_2pot, 'rays1', rays1(logical(mask_2pot), :), 'rays2', rays2(logical(mask_2pot), :), 'zeroscrewtransl', true);
                        if isfield(pose2ot, 'R')
                            if sum(mask_2pot) > sum(mask_4pst0)
                                curr_t_err_4pst0(i) = acosd(dot(gt_nt, pose2ot.t));
                            end
                        end
                    end
                    disp([sum(mask_4pst0), sum(mask_5p)]);
                    disp([gt_poses{prev_i}, gt_poses{i}]);
                    disp([gt_rel(1:3, 1:3), normc(gt_rel(1:3, 4)), gt_rel(1:3, 4), gt_poses{i}(1:3, :)]);
                    
                    if isfield(pose1, 'R'); disp([pose1.R, normc(pose1.t)]); end
                    if isfield(pose2, 'R'); disp([pose2.R, normc(pose2.t)]); end
                    disp('---');
                end
            end
            [prev_points, prev_feat, prev_im] = deal(curr_points, curr_feat, curr_im);
            prev_i = i;
        end
    end
    t_err_5p = [t_err_5p; curr_t_err_5p];
    t_err_4pst0 = [t_err_4pst0; curr_t_err_4pst0];
end

fprintf('%.1f & %.1f \\\\ \n', quantile([t_err_4pst0(:), t_err_5p(:)], [0.25, 0.5, 0.75])')
                