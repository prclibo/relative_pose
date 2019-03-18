%      1
% 
% median error
%     2.3128    2.9540    1.7819    1.7185
% 
% mean error
%     9.5477   12.4248    4.5829    3.7135
% 
% 2pot rate
%     0.3260
% 
%      2
% 
% median error
%     2.8728    3.4603    2.0810    1.9280
% 
% mean error
%     8.4015    9.9069    3.2196    2.6064
% 
% 2pot rate
%     0.2194
% 
%      3
% 
% median error
%     4.2527    4.4976    2.5972    1.7854
% 
% mean error
%    12.4082   12.3074    5.3436    2.6569
% 
% 2pot rate
%     0.0950
% 
%      4
% 
% median error
%     1.4108    1.5679    1.2196    1.6200
% 
% mean error
%     1.8296    2.2879    1.6720    2.0938
% 
% 2pot rate
%     0.4978
% 
%      5
% 
% sbisect: roots too close together
% sbisect: overflow min -2.041259 max -2.041259 diff 4.440892e-16                      nroot 2 n1 0 n2 2
% median error
%     3.0498    4.0060    2.0519    1.9873
% 
% mean error
%    11.0057   13.3671    3.5579    3.2304
% 
% 2pot rate
%     0.2292
% 
%      6
% 
% median error
%     1.6580    2.1397    1.4193    1.7902
% 
% mean error
%     5.9604    7.7781    2.2376    2.5890
% 
% 2pot rate
%     0.4042
% 
%      7
% 
% median error
%     3.9388    5.3567    2.6536    2.6651
% 
% mean error
%    16.1950   19.7457    5.7009    4.6196
% 
% 2pot rate
%     0.2292
% 
%      8
% 
% median error
%     3.0863    3.7919    2.2239    2.2500
% 
% mean error
%    12.1923   14.9149    4.7346    3.6840
% 
% 2pot rate
%     0.2339
% 
%      9
% 
% median error
%     3.4945    4.4246    2.0163    1.9673
% 
% mean error
%     9.8266   11.5348    3.6851    2.7340
% 
% 2pot rate
%     0.1676
% 
%     10
% 
% median error
%     3.6420    5.0419    2.1573    1.9340
% 
% mean error
%    12.9539   14.2233    4.0292    2.7299
% 
% 2pot rate
%     0.1876



addpath('~/workspace/relative_pose/build/matlab');
addpath('~/workspace/relative_pose/build/Matlab');
addpath('../utils');

% clear, clc
%
% seq_i = 4; % 6
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
screw_transl = nan(size(im_indices));
rotate_angle = nan(size(im_indices));
degen_2pot = false(numel(im_indices), 2);

min_move = 1;
gt_move = 0;

for i = 1:numel(im_indices)
%     fprintf('%6d / %d\n', i, numel(im_indices));
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
        
        gt_rel = normc(gt_rel);
        aaxis = real(vrrotmat2vec(gt_rel(1:3, 1:3)));
        rotate_angle(i) = aaxis(end);
        screw_transl(i) = normr(aaxis(1:3)) * normc(gt_rel(1:3, 4));
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
            
%             r4vec = vrrotmat2vec(vo_rel(1:3, 1:3));
            tic, [E_3prast0, mask_3prast0] = estimateRelativePose_PC3PRAST0_T2D(...
                aaxis(end), rays1, rays2, 0.99, thresh);
            time_3prast0(i) = toc();
            if ~isempty(E_3prast0)
                pose3 = recoverRelativePose(E_3prast0, 'rays1', rays1(logical(mask_3prast0), :), 'rays2', rays2(logical(mask_3prast0), :), 'zeroscrewtransl', true);
                if isfield(pose3, 'R')
                    t_err_3prast0(i) = acosd(dot(gt_nt, pose3.t));
                    rel_3prast0{i} = [pose3.R, pose3.t * gt_move; 0, 0, 0, 1];
                end
            end
            
            tic, [E_4pra, mask_4pra] = estimateRelativePose_PC4PRA(...
                aaxis(end), rays1, rays2, 0.99, thresh);
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
                        degen_2pot(i, 1) = true;
                    end
                    if sum(mask_2pot) > sum(mask_3prast0)
                        t_err_3prast0(i) = acosd(dot(gt_nt, pose2ot.t));
                        rel_3prast0{i} = [pose2ot.R, pose2ot.t * gt_move; 0, 0, 0, 1];
                        degen_2pot(i, 2) = true;
                    end
                end
            end
            
%             disp([sum(mask_4pst0), sum(mask_2pot), sum(mask_5p)]);
%             if isfield(pose1, 'R'); disp([pose1.R, normc(pose1.t)]); end
%             if isfield(pose2, 'R'); disp([pose2.R, normc(pose2.t)]); end
%             disp([pose2ot.R, normc(pose2ot.t)]);
%             disp([pose3.R, normc(pose3.t)]);
%             disp([pose4.R, normc(pose4.t)]);
%             disp([gt_rel(1:3, 1:3), normc(gt_rel(1:3, 4)), gt_rel(1:3, 4)]);
%             disp('---');

%             return;
%             if i == 11
%                 return; 
%             end
            
        end
        [prev_points, prev_feat, prev_im] = deal(curr_points, curr_feat, curr_im);
        prev_i = i;
%         if i > 300; return; end
    end
end

% save(sprintf('/tmp/ws_kitti_eq%02d_move%d', seq_i, min_move));