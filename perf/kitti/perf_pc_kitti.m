addpath('~/workspace/relative_pose/build/matlab');
addpath('~/workspace/relative_pose/build/Matlab');
addpath('../utils');

clear, clc

seq_i = 5; %1; %4;
% 6? 9? 12?
cam_dir = sprintf('~/data/kitti/data_odometry_gray/sequences/%02d', seq_i);
pose_path = sprintf('~/data/kitti/dataset/poses/%02d.txt', seq_i);

im_dir = fullfile(cam_dir, 'image_0');
calib_path = fullfile(cam_dir, 'calib.txt');
calib = loadCalib(calib_path);
calib = calib(1);

gt_pose_data = load(pose_path);
gt_poses = cell(size(gt_pose_data, 1), 1);
for i = 1:numel(gt_poses)
    gt_poses{i} = reshape(gt_pose_data(i, :), 4, 3)';
    gt_poses{i}(4, 4) = 1;
end

fx = calib.intrs(1, 1);
fy = calib.intrs(2, 2);
cx = calib.intrs(1, 3);
cy = calib.intrs(2, 3);
thresh = 2 / fx;

consec_interv = 1;
im_indices = 20:consec_interv:numel(gt_poses) - 20;

t_err_5p = nan(size(im_indices));
t_err_4pst0 = nan(size(im_indices));
t_err_4pra = nan(size(im_indices));
t_err_3prast0 = nan(size(im_indices));
rel_5p = cell(size(im_indices));
rel_4pst0 = cell(size(im_indices));
rel_4pra = cell(size(im_indices));
rel_3prast0 = cell(size(im_indices));
time_5p = nan(size(im_indices));
time_4pst0 = nan(size(im_indices));
time_4pra = nan(size(im_indices));
time_3prast0 = nan(size(im_indices));

min_move = 1;
gt_move = 0;

for i = 1:numel(im_indices)
    fprintf('%6d / %d\n', i, numel(im_indices));
    rng(3);
%     clc;
    if i > 1
        offseted = gt_poses{i};
        offseted_prev = gt_poses{prev_i};
        offseted(1:3, 4) = offseted(1:3, 4) - offseted_prev(1:3, 4);
        offseted_prev(1:3, 4) = 0;
        gt_rel = offseted_prev \ offseted;
        gt_rel = calib.extrs \ gt_rel * calib.extrs;
        gt_nt = normc(gt_rel(1:3, 4));
        gt_move = norm(gt_rel(1:3, 4));
    end
    if i == 1 || gt_move > min_move
        file_name = sprintf('%06d.png', im_indices(i));
        curr_path = fullfile(im_dir, file_name);
        curr_im = imread(curr_path);
        curr_im = histeq(curr_im, 255);
        curr_points = detectSURFFeatures(curr_im, ...
            'metricthreshold', 300); % , 'roi', im_roi);
        [curr_feat, curr_points] = extractFeatures(curr_im, curr_points);
        if i > 1% && gt_move > min_move
            pairs = matchFeatures(curr_feat, prev_feat,...
                'method', 'Exhaustive', 'maxratio', 0.6);
            if size(pairs, 1) < 50; continue; end
            
            mpoints1 = curr_points(pairs(:, 1));
            mpoints2 = prev_points(pairs(:, 2));
%             figure; showMatchedFeatures(curr_im, prev_im, mpoints1, mpoints2);
            rays1 = ([mpoints1(:).Location] - [cx, cy]) ./ [fx, fy];
            rays2 = ([mpoints2(:).Location] - [cx, cy]) ./ [fx, fy];
            rays1(:, 3) = 1; rays2(:, 3) = 1;
                        
            tic, [E_5p, mask_5p] = estimateRelativePose_PC5P_LiH(...
                rays1, rays2, 0.99, thresh);
            time_5p(i) = toc();
            if ~isempty(E_5p)
                pose1 = recoverRelativePose(E_5p, 'rays1', rays1(logical(mask_5p), :), 'rays2', rays2(logical(mask_5p), :));
                if (isfield(pose1, 'R'))
                    t_err_5p(i) = acosd(dot(gt_nt, pose1.t));
                    rel_5p{i} = [pose1.R, pose1.t * gt_move; 0, 0, 0, 1];
                end
            end
            
            r4vec = vrrotmat2vec(gt_rel(1:3, 1:3));
%             r4vec = vrrotmat2vec(vo_rel(1:3, 1:3));
            tic, [E_3prast0, mask_3prast0] = estimateRelativePose_PC3PRAST0_T2D(...
                r4vec(end), rays1, rays2, 0.99, thresh);
            time_3prast0(i) = toc();
            if ~isempty(E_3prast0)
                pose3 = recoverRelativePose(E_3prast0, 'rays1', rays1(logical(mask_3prast0), :), 'rays2', rays2(logical(mask_3prast0), :), 'zeroscrewtransl', true);
                if isfield(pose3, 'R')
                    t_err_3prast0(i) = acosd(dot(gt_nt, pose3.t));
                    rel_3prast0{i} = [pose3.R, pose3.t * gt_move; 0, 0, 0, 1];
                end
            end
            
            tic, [E_4pra, mask_4pra] = estimateRelativePose_PC4PRA(...
                r4vec(end), rays1, rays2, 0.99, thresh);
            time_4pra(i) = toc();
            if ~isempty(E_4pra)
                pose4 = recoverRelativePose(E_4pra, 'rays1', rays1(logical(mask_4pra), :), 'rays2', rays2(logical(mask_4pra), :), 'zeroscrewtransl', false);
                if isfield(pose4, 'R')
                    t_err_4pra(i) = acosd(dot(gt_nt, pose4.t));
                    rel_4pra{i} = [pose4.R, pose4.t * gt_move; 0, 0, 0, 1];
                end
            end
            
            tic, [E_4pst0, mask_4pst0] = estimateRelativePose_PC4PST0_NullE(...
                rays1, rays2, 0.99, thresh);
            time_4pst0(i) = toc();
            if ~isempty(E_4pst0)
                pose2 = recoverRelativePose(E_4pst0, 'rays1', rays1(logical(mask_4pst0), :), 'rays2', rays2(logical(mask_4pst0), :), 'zeroscrewtransl', true);
                if isfield(pose2, 'R')
                    t_err_4pst0(i) = acosd(dot(gt_nt, pose2.t));
                    rel_4pst0{i} = [pose2.R, pose2.t * gt_move; 0, 0, 0, 1];
                end
            end
            [E_2pot, mask_2pot] = estimateRelativePose_PC2POT(...
                rays1, rays2, 0.99, thresh);
            if ~isempty(E_2pot)
                pose2ot = recoverRelativePose(E_2pot, 'rays1', rays1(logical(mask_2pot), :), 'rays2', rays2(logical(mask_2pot), :), 'zeroscrewtransl', true);
                if isfield(pose2ot, 'R')
                    if sum(mask_2pot) > sum(mask_4pst0)
                        t_err_4pst0(i) = acosd(dot(gt_nt, pose2ot.t));
                        rel_4pst0{i} = [pose2ot.R, pose2ot.t * gt_move; 0, 0, 0, 1];
                    end
                    if sum(mask_2pot) > sum(mask_3prast0)
                        t_err_3prast0(i) = acosd(dot(gt_nt, pose2ot.t));
                        rel_3prast0{i} = [pose2ot.R, pose2ot.t * gt_move; 0, 0, 0, 1];
                    end
                end
            end
            
            disp([sum(mask_4pst0), sum(mask_2pot), sum(mask_5p)]);
            disp([pose1.R, normc(pose1.t)]);
            disp([pose2.R, normc(pose2.t)]);
%             disp([pose2ot.R, normc(pose2ot.t)]);
%             disp([pose3.R, normc(pose3.t)]);
%             disp([pose4.R, normc(pose4.t)]);
            disp([gt_rel(1:3, 1:3), normc(gt_rel(1:3, 4)), gt_rel(1:3, 4)]);
            disp('---');

%             return;
%             if i == 11
%                 return; 
%             end
            
        end
        [prev_points, prev_feat, prev_im] = deal(curr_points, curr_feat, curr_im);
        prev_i = i;
        if i > 200; return; end
    end
end