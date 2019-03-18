clear, clc

rng(3);

sdk_dir = '~/workspace/robotcar-dataset-sdk';
% data_dir = '~/data/robotcar/2014-05-14-13-46-12';
data_dir = '~/data/robotcar/2014-06-26-08-53-56';
model_dir = fullfile(sdk_dir, 'models');
extrinsic_dir = fullfile(sdk_dir, 'extrinsics');
ins_file = fullfile(data_dir, 'gps/ins.csv');
vo_file = fullfile(data_dir, 'vo/vo.csv');

addpath(fullfile(sdk_dir, 'matlab'));
addpath(fullfile(fileparts(mfilename('fullpath')), '/../../build/matlab/'));

% error('consec = 1 with choice 2 can outperform');
choice = 2;
if choice == 1
    image_dir = fullfile(data_dir, 'stereo/centre');
    timestamp_file = fullfile(data_dir, 'stereo.timestamps');
    ins_extrs = SE3MatrixFromComponents(dlmread([extrinsic_dir, '/ins.txt']));
    cam_extrs = SE3MatrixFromComponents(dlmread([extrinsic_dir, '/stereo.txt']));
    [fx, fy, cx, cy, G_camera_image, LUT] = ReadCameraModel(image_dir, model_dir);
    cam_ins = ins_extrs \ cam_extrs * G_camera_image;
    lower_clip = 830;
elseif choice == 2
    image_dir = fullfile(data_dir, 'mono_rear');
    timestamp_file = fullfile(data_dir, 'mono_rear.timestamps');
    ins_extrs = SE3MatrixFromComponents(dlmread([extrinsic_dir, '/ins.txt']));
    cam_extrs = SE3MatrixFromComponents(dlmread([extrinsic_dir, '/mono_rear.txt']));
    [fx, fy, cx, cy, G_camera_image, LUT] = ReadCameraModel(image_dir, model_dir);
    cam_ins = ins_extrs \ cam_extrs * G_camera_image;
    lower_clip = 950;
end

consec_interv = 1;
timestamps = load(timestamp_file);
timestamps = timestamps(1:consec_interv:end, 1)';

thresh = 4 / fx;

scale = 0.5;
fx = fx * scale;
fy = fy * scale;
cx = cx * scale;
cy = cy * scale;
thresh = thresh * scale;
lower_clip = lower_clip * scale;

% timestamps = timestamps(197:200);
timestamps = timestamps(50:end - 50);

ins_poses = InterpolatePoses(ins_file, timestamps, timestamps(1));
vo_poses = RelativeToAbsolutePoses(vo_file, timestamps, timestamps(1));
        
t_err_5p = nan(size(timestamps));
t_err_4pst0 = nan(size(timestamps));
t_err_4pra = nan(size(timestamps));
t_err_3prast0 = nan(size(timestamps));
rel_5p = cell(size(timestamps));
rel_4pst0 = cell(size(timestamps));
rel_4pra = cell(size(timestamps));
rel_3prast0 = cell(size(timestamps));
time_5p = nan(size(timestamps));
time_4pst0 = nan(size(timestamps));
time_4pra = nan(size(timestamps));
time_3prast0 = nan(size(timestamps));

min_move = 2;
vo_move = 0;
for i = 1:numel(timestamps)
    fprintf('%6d / %d\n', i, numel(timestamps));
    rng(3);
%     clc;
    
    if i > 1
        offseted = vo_poses{i};
        offseted_prev = vo_poses{prev_i};
        offseted(1:3, 4) = offseted(1:3, 4) - offseted_prev(1:3, 4);
        offseted_prev(1:3, 4) = 0;
        vo_rel = offseted_prev \ offseted;
        vo_rel = cam_ins \ vo_rel * cam_ins;
        ins_rel = ins_poses{prev_i} \ ins_poses{i};
        vo_nt = normc(vo_rel(1:3, 4));
        vo_move = norm(vo_rel(1:3, 4));
    end
    if i == 1 || vo_move > min_move
        curr_im = rgb2gray(LoadImage(image_dir, timestamps(i), LUT));
        curr_im = imresize(curr_im, scale);
        curr_im = curr_im(1:lower_clip, :);
        temp = curr_im;
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
            
            if isempty(rays1) || isempty(rays2); continue; end
            
            tic, [E_5p, mask_5p] = estimateRelativePose_PC5P_LiH(...
                rays1, rays2, 0.999, thresh);
            time_5p(i) = toc();
            if ~isempty(E_5p)
                pose1 = recoverRelativePose(E_5p, 'rays1', rays1(logical(mask_5p), :), 'rays2', rays2(logical(mask_5p), :));
                if (isfield(pose1, 'R'))
                    t_err_5p(i) = acosd(dot(vo_nt, pose1.t));
                    rel_5p{i} = [pose1.R, pose1.t * vo_move; 0, 0, 0, 1];
                end
            end
            
            tic, [E_4pst0, mask_4pst0] = estimateRelativePose_PC4PST0_NullE(...
                rays1, rays2, 0.999, thresh);
            time_4pst0(i) = toc();
            if ~isempty(E_4pst0)
                pose2 = recoverRelativePose(E_4pst0, 'rays1', rays1(logical(mask_4pst0), :), 'rays2', rays2(logical(mask_4pst0), :), 'zeroscrewtransl', true);
                if isfield(pose2, 'R')
                    t_err_4pst0(i) = acosd(dot(vo_nt, pose2.t));
                    rel_4pst0{i} = [pose2.R, pose2.t * vo_move; 0, 0, 0, 1];
                end
            end
            
            r4vec = vrrotmat2vec(ins_rel(1:3, 1:3));
%             r4vec = vrrotmat2vec(vo_rel(1:3, 1:3));
            tic, [E_3prast0, mask_3prast0] = estimateRelativePose_PC3PRAST0_T2D(...
                r4vec(end), rays1, rays2, 0.999, thresh);
            time_3prast0(i) = toc();
            if ~isempty(E_3prast0)
                pose3 = recoverRelativePose(E_3prast0, 'rays1', rays1(logical(mask_3prast0), :), 'rays2', rays2(logical(mask_3prast0), :), 'zeroscrewtransl', true);
                if isfield(pose3, 'R')
                    t_err_3prast0(i) = acosd(dot(vo_nt, pose3.t));
                    rel_3prast0{i} = [pose3.R, pose3.t * vo_move; 0, 0, 0, 1];
                end
            end
            
            tic, [E_4pra, mask_4pra] = estimateRelativePose_PC4PRA(...
                r4vec(end), rays1, rays2, 0.999, thresh);
            time_4pra(i) = toc();
            if ~isempty(E_4pra)
                pose4 = recoverRelativePose(E_4pra, 'rays1', rays1(logical(mask_4pra), :), 'rays2', rays2(logical(mask_4pra), :), 'zeroscrewtransl', true);
                if isfield(pose4, 'R')
                    t_err_4pra(i) = acosd(dot(vo_nt, pose4.t));
                    rel_4pra{i} = [pose4.R, pose4.t * vo_move; 0, 0, 0, 1];
                end
            end
            
            disp([sum(mask_4pst0), sum(mask_5p)]);
            disp([pose1.R, normc(pose1.t)]);
            disp([pose2.R, normc(pose2.t)]);
            disp([pose3.R, normc(pose3.t)]);
            disp([pose4.R, normc(pose4.t)]);
            disp([vo_rel(1:3, 1:3), normc(vo_rel(1:3, 4)), vo_rel(1:3, 4)]);
            disp('---');
        end
        [prev_points, prev_feat, prev_im] = deal(curr_points, curr_feat, curr_im);
        prev_i = i;
        if i > 200; return; end
    end 
end
fprintf('\n');