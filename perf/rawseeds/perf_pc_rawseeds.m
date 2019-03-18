
% Bicocca_2009-02-26a odometry
%     7.3173   10.5614    7.0668    9.1367
%    14.6692   18.7773   12.5030   15.1218
%    27.4657   33.7739   21.7977   23.6057
% % Bicocca_2009-02-26a imu
%     4.5205    6.7192    7.0668    9.1892
%     7.4796   11.0883   12.5255   15.2482
%    11.5709   15.7383   21.7977   23.6116

clear;
rng(3);
rot_type = 'imu';

addpath('~/data/umich_ford/Code/MATLAB/create_ijrr_utils/')
addpath(fullfile(fileparts(mfilename('fullpath')), '/../../build/matlab/'));
addpath(fullfile(fileparts(mfilename('fullpath')), '/../utils'));

data_name = 'Bicocca_2009-02-26a';
data_dir = fullfile('~/data/rawseeds', data_name);

im_dir = fullfile(data_dir, 'FRONTAL');
gt_path = fullfile(data_dir, sprintf('%s-GROUNDTRUTH.csv', data_name));
imu_path = fullfile(data_dir, sprintf('%s-IMU_STRETCHED.csv', data_name));
odo_path = fullfile(data_dir, sprintf('%s-ODOMETRY_XYT.csv', data_name));
im_list_path = fullfile(data_dir, sprintf('%s-LISTS/%s-FRONTAL.lst', data_name, data_name));
stamp_path = fullfile(data_dir, sprintf('%s-LISTS/%s-FRONTAL.csv', data_name, data_name));
calib_path = fullfile(data_dir, '../Calibration_04-Results/Calibration_04-Intrinsics_FRONTAL.mat');

fid_im_list = fopen(im_list_path, 'r');
im_list = textscan(fid_im_list, '%s');
im_list = im_list{1};
fclose(fid_im_list);

cam_intrs = load(calib_path);
im_stamps = load(stamp_path);
[gt_poses, im_indices] = loadRawSeedsGTPoses(gt_path, im_stamps);
im_stamps = im_stamps(im_indices);
im_list = im_list(im_indices);

if strcmp(rot_type, 'odo')
    [odo_poses, im_indices] = loadRawSeedsOdometry(odo_path, im_stamps);
elseif strcmp(rot_type, 'imu')
    [imu_rotations, im_indices] = loadRawSeedsIMURotations(imu_path, im_stamps);
else; error('Wrong rot type');
end
gt_poses = gt_poses(im_indices);
im_stamps = im_stamps(im_indices);
im_list = im_list(im_indices);

cam_extrs = [ 0,  0, 1, -0.171;
             -1,  0, 0, -0.115;
              0, -1, 0,  1.068;
              0,  0, 0,  1];
odo_extrs = [ 0, 1, 0, 0;
             -1, 0, 0, 0;
              0, 0, 1, 0;
              0, 0, 0, 1];

fx = cam_intrs.fc_frontal(1);
fy = cam_intrs.fc_frontal(2);
cx = cam_intrs.cc_frontal(1);
cy = cam_intrs.cc_frontal(2);

cam_param = cameraParameters(...
    'IntrinsicMatrix', [fx, 0, 0; 0, fy, 0; cx, cy, 1], ...
    'RadialDistortion', cam_intrs.kc_frontal(1:2));

thresh = 1 / fx;

t_err_5p = nan(size(im_stamps));
t_err_4pst0 = nan(size(im_stamps));
t_err_4pra = nan(size(im_stamps));
t_err_3prast0 = nan(size(im_stamps));
rel_5p = cell(size(im_stamps));
rel_4pst0 = cell(size(im_stamps));
rel_4pra = cell(size(im_stamps));
rel_3prast0 = cell(size(im_stamps));
time_5p = nan(size(im_stamps));
time_4pst0 = nan(size(im_stamps));
time_4pra = nan(size(im_stamps));
time_3prast0 = nan(size(im_stamps));

min_move = 0.1;
gt_move = 0;
for i = 1:numel(im_stamps)
    fprintf('%d / %d\n', i, numel(im_stamps))
    rng(3);
%     clc
    if i > 1
        gt_rel = relativePose(gt_poses{i}, gt_poses{prev_i});
        gt_rel = cam_extrs \ gt_rel * cam_extrs;
        if strcmp(rot_type, 'odo')
            odo_rel = relativePose(odo_poses{i}, odo_poses{prev_i});
            %         odo_rel = odo_extrs * odo_rel * inv(odo_extrs);
            odo_rel = cam_extrs \ odo_rel * cam_extrs;
            rot_rel = odo_rel(1:3, 1:3);
        elseif strcmp(rot_type, 'imu')
            imu_rel = imu_rotations{prev_i} \ imu_rotations{i};
            rot_rel = imu_rel;
        else; error('Wrong type');
        end
        
        gt_nt = normc(gt_rel(1:3, 4));
        gt_move = norm(gt_rel(1:3, 4));
        gt_E = skew(gt_nt) * gt_rel(1:3, 1:3);
        gt_E = gt_E / norm(gt_E(:));
    end
    if i == 1 || gt_move > min_move
        curr_path = fullfile(im_dir, im_list{i});
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
                
                r4vec = vrrotmat2vec(rot_rel(1:3, 1:3));
                %             r4vec = vrrotmat2vec(vo_rel(1:3, 1:3));
                tic, [E_3prast0, mask_3prast0] = estimateRelativePose_PC3PRAST0_T2D(...
                    r4vec(end), rays1, rays2, 0.999, thresh);
                time_3prast0(i) = toc();
                if ~isempty(E_3prast0)
                    pose3 = recoverRelativePose(E_3prast0, 'rays1', rays1(logical(mask_3prast0), :), 'rays2', rays2(logical(mask_3prast0), :), 'zeroscrewtransl', true);
                    if isfield(pose3, 'R')
                        t_err_3prast0(i) = acosd(dot(gt_nt, pose3.t));
                        rel_3prast0{i} = [pose3.R, pose3.t * gt_move; 0, 0, 0, 1];
                    end
                end
                
                tic, [E_4pra, mask_4pra] = estimateRelativePose_PC4PRA(...
                    r4vec(end), rays1, rays2, 0.999, thresh);
                time_4pra(i) = toc();
                if ~isempty(E_4pra)
                    pose4 = recoverRelativePose(E_4pra, 'rays1', rays1(logical(mask_4pra), :), 'rays2', rays2(logical(mask_4pra), :), 'zeroscrewtransl', true);
                    if isfield(pose4, 'R')
                        t_err_4pra(i) = acosd(dot(gt_nt, pose4.t));
                        rel_4pra{i} = [pose4.R, pose4.t * gt_move; 0, 0, 0, 1];
                    end
                end
                
                [E_2pot, mask_2pot] = estimateRelativePose_PC2POT(...
                    rays1, rays2, 0.999, thresh);
%                 if ~isempty(E_2pot)
%                     pose2ot = recoverRelativePose(E_2pot, 'rays1', rays1(logical(mask_2pot), :), 'rays2', rays2(logical(mask_2pot), :), 'zeroscrewtransl', true);
%                     if isfield(pose2ot, 'R')
%                         if sum(mask_2pot) > sum(mask_4pst0)
%                             t_err_4pst0(i) = acosd(dot(gt_nt, pose2ot.t));
%                             rel_4pst0{i} = [pose2ot.R, pose2ot.t * gt_move; 0, 0, 0, 1];
%                         end
%                         if sum(mask_2pot) > sum(mask_3prast0)
%                             t_err_3prast0(i) = acosd(dot(gt_nt, pose2ot.t));
%                             rel_3prast0{i} = [pose2ot.R, pose2ot.t * gt_move; 0, 0, 0, 1];
%                         end
%                     end
%                 end

                
                disp([sum(mask_4pst0), sum(mask_5p)]);
                disp([gt_rel(1:3, 1:3), normc(gt_rel(1:3, 4)), gt_rel(1:3, 4), gt_E]);

                if isfield(pose1, 'R'); disp([pose1.R, normc(pose1.t)]); end
                if isfield(pose2, 'R'); disp([pose2.R, normc(pose2.t)]); end
                %             disp([pose2ot.R, normc(pose2ot.t)]);
                %             disp([pose3.R, normc(pose3.t)]);
                %             disp([pose4.R, normc(pose4.t)]);
                disp('---');

            end
        end
        [prev_points, prev_feat, prev_im] = deal(curr_points, curr_feat, curr_im);
        prev_i = i;
%         if i > 3000; return; end
    end
end