clear;

rng(3);

sdk_dir = '~/workspace/robotcar-dataset-sdk';
data_dir = '~/data/robotcar/temp/2014-05-14-13-46-12';
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
thresh = 8 / fx;
ins_extrs = SE3MatrixFromComponents(dlmread([extrinsic_dir, '/ins.txt']));
stereo_extrs = SE3MatrixFromComponents(dlmread([extrinsic_dir, '/stereo.txt']));
camera_ins = ins_extrs \ stereo_extrs * G_camera_image;

timestamps = timestamps(197:200);
otime = timestamps(1);
ins_poses = InterpolatePoses(ins_file, timestamps, otime);
vo_poses = RelativeToAbsolutePoses(vo_file, timestamps, otime);

t_err_5p = nan(numel(timestamps), 1);
t_err_4pst0 = nan(numel(timestamps), 1);
% fprintf('Processing          ');
for i = 1:numel(timestamps)
%     fprintf('\b\b\b\b\b\b%6d', i);
    rng(3);

    curr_im = rgb2gray(LoadImage(image_dir, timestamps(i), LUT));
    curr_points = detectSURFFeatures(curr_im);
    [curr_feat, curr_points] = extractFeatures(curr_im, curr_points);
    
    if i > 1
        pairs = matchFeatures(curr_feat, prev_feat); %ZZG"JGJFFFFFF, 'method', 'approximate') ;
        mpoints1 = curr_points(pairs(:, 1));
        mpoints2 = prev_points(pairs(:, 2));
        
        dist = sqrt(sum(([mpoints1(:).Location] - [mpoints2(:).Location]).^2, 2));
        mpoints1 = mpoints1(dist > 10);
        mpoints2 = mpoints2(dist > 10);
        
        rays1 = ([mpoints1(:).Location] - [cx, cy]) ./ [fx, fy];
        rays2 = ([mpoints2(:).Location] - [cx, cy]) ./ [fx, fy];
        rays1(:, 3) = 1; rays2(:, 3) = 1;
        dist = sqrt(sum(([mpoints1(:).Location] - [mpoints2(:).Location]).^2, 2));
                
        vo_relative = camera_ins \ (vo_poses{i - 1} \ vo_poses{i}) * camera_ins;
        ins_relative = camera_ins \ (ins_poses{i - 1} \ ins_poses{i}) * camera_ins;
        vo_t = normc(vo_relative(1:3, 4));
        ins_t = normc(ins_relative(1:3, 4));
        
        ins_pose.R = ins_relative(1:3, 1:3);
        ins_pose.t = ins_relative(1:3, 4);

        [E_5p, mask_5p] = estimateRelativePose_PC5P_LiH(rays1, rays2, 0.999, thresh);
        fprintf('mask %f\n', sum(mask_5p));
        if ~isempty(E_5p)
            pose1 = recoverRelativePose(E_5p, 'rays1', rays1(logical(mask_5p), :), 'rays2', rays2(logical(mask_5p), :));
            t_err_5p(i) = acosd(dot(vo_t, pose1.t));
        end

        [E_4pst0, mask_4pst0] = estimateRelativePose_PC4PST0_NullE(rays1, rays2, 0.999, thresh);
        if ~isempty(E_4pst0)
            pose2 = recoverRelativePose(E_4pst0, 'rays1', rays1(logical(mask_4pst0), :), 'rays2', rays2(logical(mask_4pst0), :), 'zeroscrewtransl', true);
            t_err_4pst0(i) = acosd(dot(vo_t, pose2.t));
        end

        [E_2pot, mask_2pot] = estimateRelativePose_PC2POT(rays1, rays2, 0.999, thresh);
        if ~isempty(E_2pot)
            pose3 = recoverRelativePose(E_2pot, 'rays1', rays1(logical(mask_2pot), :), 'rays2', rays2(logical(mask_2pot), :), 'zeroscrewtransl', true);
            if sum(mask_2pot) >= sum(mask_4pst0) * 0.99
                t_err_4pst0(i) = acosd(dot(vo_t, pose3.t));
            end
        end
        
        t_vo_ins(i) = acosd(dot(vo_t, ins_t));
%         figure; showMatchedFeatures(curr_im, prev_im, mpoints1, mpoints2);
        
        fprintf('time: %f\n', timestamps(i));
        fprintf('vo_t: %f, %f, %f\n', vo_t);
        fprintf('ins_t: %f, %f, %f\n', ins_t);
        fprintf('pose1.t: %f, %f, %f\n', pose1.t);
        fprintf('pose2.t: %f, %f, %f\n', pose2.t);
        fprintf('pose3.t: %f, %f, %f\n', pose3.t);

%         return;
    end 
    [prev_points, prev_feat, prev_im] = deal(curr_points, curr_feat, curr_im);
end
fprintf('\n');