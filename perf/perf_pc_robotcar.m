clear, clc

rng(3);

sdk_dir = '~/workspace/robotcar-dataset-sdk';
data_dir = '~/data/robotcar/2014-05-14-13-46-12';
data_dir = '~/data/robotcar/2014-06-26-08-53-56';
image_dir = fullfile(data_dir, 'stereo/centre');
timestamp_file = fullfile(data_dir, 'stereo.timestamps');
model_dir = fullfile(sdk_dir, 'models');
extrinsic_dir = fullfile(sdk_dir, 'extrinsics');
ins_file = fullfile(data_dir, 'gps/ins.csv');
vo_file = fullfile(data_dir, 'vo/vo.csv');

addpath(fullfile(sdk_dir, 'matlab'));
addpath(fullfile(fileparts(mfilename('fullpath')), '/../build/matlab/'));

consec_interv = 5;
timestamps = load(timestamp_file);
timestamps = timestamps(1:consec_interv:end, 1)';
[fx, fy, cx, cy, G_camera_image, LUT] = ReadCameraModel(image_dir, model_dir);
thresh = 6 / fx;
ins_extrs = SE3MatrixFromComponents(dlmread([extrinsic_dir, '/ins.txt']));
stereo_extrs = SE3MatrixFromComponents(dlmread([extrinsic_dir, '/stereo.txt']));
camera_ins = ins_extrs \ stereo_extrs * G_camera_image;
lower_clip = 830;

scale = 0.5;
fx = fx * scale;
fy = fy * scale;
cx = cx * scale;
cy = cy * scale;
thresh = thresh * scale;
lower_clip = lower_clip * scale;

% timestamps = timestamps(197:200);
timestamps = timestamps(52 + 75:250);

t_err_5p = nan(numel(timestamps), 1);
t_err_4pst0 = nan(numel(timestamps), 1);
rel_5p = cell(numel(timestamps), 1);
rel_4pst0 = cell(numel(timestamps), 1);

% fprintf('Processing          ');
for i = 1:numel(timestamps)
    fprintf('%6d / %d\n', i, numel(timestamps));
    rng(3);

    curr_im = rgb2gray(LoadImage(image_dir, timestamps(i), LUT));
    curr_im = imresize(curr_im, scale);
    curr_im = curr_im(1:lower_clip, :);
    curr_points = detectSURFFeatures(curr_im);
    [curr_feat, curr_points] = extractFeatures(curr_im, curr_points);
    
    if i > 1
        ins_pose = InterpolatePoses(ins_file, timestamps(i), timestamps(i - 1));
        vo_pose = RelativeToAbsolutePoses(vo_file, timestamps(i), timestamps(i - 1));
        
        pairs = matchFeatures(curr_feat, prev_feat,...
            'method', 'Exhaustive', 'maxratio', 0.6) ;
        mpoints1 = curr_points(pairs(:, 1));
        mpoints2 = prev_points(pairs(:, 2));
        
%         dist = sqrt(sum(([mpoints1(:).Location] - [mpoints2(:).Location]).^2, 2));
%         mpoints1 = mpoints1(dist > 10);
%         mpoints2 = mpoints2(dist > 10);
        
        rays1 = ([mpoints1(:).Location] - [cx, cy]) ./ [fx, fy];
        rays2 = ([mpoints2(:).Location] - [cx, cy]) ./ [fx, fy];
        rays1(:, 3) = 1; rays2(:, 3) = 1;
        dist = sqrt(sum(([mpoints1(:).Location] - [mpoints2(:).Location]).^2, 2));
                
        vo_relative = camera_ins \ vo_pose{1} * camera_ins;
        ins_relative = camera_ins \ ins_pose{1} * camera_ins;
        vo_nt = normc(vo_relative(1:3, 4));
        ins_nt = normc(ins_relative(1:3, 4));

        [E_5p, mask_5p] = estimateRelativePose_PC5P_LiH(rays1, rays2, 0.999, thresh);
%         fprintf('mask %f\n', sum(mask_5p));
        if ~isempty(E_5p)
            pose1 = recoverRelativePose(E_5p, 'rays1', rays1(logical(mask_5p), :), 'rays2', rays2(logical(mask_5p), :));
            if (isfield(pose1, 'R'))
                t_err_5p(i) = acosd(dot(vo_nt, pose1.t));
                rel_5p{i} = [pose1.R, pose1.t; 0, 0, 0, 1];
            end
        end

        [E_4pst0, mask_4pst0] = estimateRelativePose_PC4PST0_NullE(rays1, rays2, 0.999, thresh);
        if ~isempty(E_4pst0)
            pose2 = recoverRelativePose(E_4pst0, 'rays1', rays1(logical(mask_4pst0), :), 'rays2', rays2(logical(mask_4pst0), :), 'zeroscrewtransl', true);
            if isfield(pose2, 'R')
                t_err_4pst0(i) = acosd(dot(vo_nt, pose2.t));
                rel_4pst0{i} = [pose2.R, pose2.t; 0, 0, 0, 1];
            end
        end

        [E_2pot, mask_2pot] = estimateRelativePose_PC2POT(rays1, rays2, 0.999, thresh);
        if ~isempty(E_2pot)
            pose3 = recoverRelativePose(E_2pot, 'rays1', rays1(logical(mask_2pot), :), 'rays2', rays2(logical(mask_2pot), :), 'zeroscrewtransl', true);
            if sum(mask_2pot) >= sum(mask_4pst0) * 0.99 && isfield(pose3, 'R')
                t_err_4pst0(i) = acosd(dot(vo_nt, pose3.t));
                rel_4pst0{i} = [pose3.R, pose3.t; 0, 0, 0, 1];
            end
        end
%         figure; showMatchedFeatures(curr_im, prev_im, mpoints1, mpoints2);

        tnorm = norm(vo_relative(1:3, 4));
        if ~isempty(rel_4pst0{i})
            rel_4pst0{i}(1:3, 4) = rel_4pst0{i}(1:3, 4) * tnorm;
        end
        if ~isempty(rel_5p{i})
            rel_5p{i}(1:3, 4) = rel_5p{i}(1:3, 4) * tnorm;
        end
        
        axis = vrrotmat2vec(ins_relative(1:3, 1:3));
        ins_st(i) = dot(axis(1:3), ins_nt);
        
        fprintf('%d time: %f\n', i - 2, timestamps(i));
        fprintf('%f %f\n', t_err_4pst0(i), t_err_5p(i));
        fprintf('%f\n', ins_st(i));

        disp(vo_relative);
        disp([pose2.R, pose2.t]);

%         return;
    end 
    [prev_points, prev_feat, prev_im] = deal(curr_points, curr_feat, curr_im);
end
fprintf('\n');