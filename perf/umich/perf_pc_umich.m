rng(3);

sdk_dir = '/Users/li/data/umich_ford/Code/MATLAB';
db_path = fullfile(sdk_dir, 'db.xml');
data_dir = '/Users/li/data/umich_ford/IJRR-Dataset-1';
im_dir = fullfile(data_dir, 'IMAGES/Cam0');
param_path = fullfile(data_dir, 'PARAM.mat');
pose_path = fullfile(data_dir, 'Pose-Applanix.log');
stamp_path = fullfile(data_dir, 'Timestamp.log');

addpath(fullfile(sdk_dir, 'create_ijrr_utils'));
addpath(fullfile(sdk_dir, 'create_ijrr_utils/fcns_coordxform'));
addpath(fullfile(sdk_dir, 'create_ijrr_utils/Multiprod-3'));
addpath(fullfile(sdk_dir, 'xmltools'));
addpath(fullfile(sdk_dir, 'xml_toolbox'));
addpath(sdk_dir);
addpath('../build/matlab')

if ~exist('db', 'var'); db = hdl_newloaddb(db_path); end
if ~exist('lcalib', 'var'); lcalib = hdl_lasergeom(db); end
if ~exist('ccalib', 'var'); ccalib = load(param_path); end
if ~exist('gt_poses', 'var'); gt_poses = load_pose_applanix(pose_path); end
if ~exist('im_list', 'var')
    fid_stamp = fopen(stamp_path, 'r');
    im_list = textscan(fid_stamp, '%f %f %f %f', 'HeaderLines', 1);
end

consec_interv = 1;
start = 500;
im_stamps = im_list{2}(start:consec_interv:end);
im_indices = im_list{1}(start:consec_interv:end);
im_stamps = im_stamps(180:500);
im_indices = im_indices(180:500);

if ~exist('cam_poses', 'var')
    cam_poses = loadCameraPoses(im_stamps, gt_poses, lcalib, ccalib.PARAM(1));
end

fx = ccalib.PARAM(1).K(1, 1);
fy = ccalib.PARAM(1).K(2, 2);
cx = ccalib.PARAM(1).K(1, 3);
cy = ccalib.PARAM(1).K(2, 3);
thresh = 2 / fx;
im_size = [1232,1616];
im_roi = [[600, 200], im_size([2, 1]) - [600, 200] * 2];

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

min_move = 0.5;
gt_move = 0;
for i = 1:numel(im_stamps)
    fprintf('%d / %d\n', i, numel(im_stamps))
    rng(3);
%     clc
    if i > 1
        offseted = cam_poses(:, :, i);
        offseted(1:3, 4) = offseted(1:3, 4) - cam_poses(1:3, 4, prev_i);
        offseted_prev = cam_poses(:, :, prev_i);
        offseted_prev(1:3, 4) = 0;
        gt_rel = offseted_prev \ offseted;
        gt_nt = normc(gt_rel(1:3, 4));
        gt_move = norm(gt_rel(1:3, 4));
        gt_E = skew(gt_nt) * gt_rel(1:3, 1:3);
        gt_E = gt_E / norm(gt_E(:));
    end
    if i == 1 || gt_move > min_move
        file_name = sprintf('image%04d.ppm', im_indices(i));
        curr_path = fullfile(im_dir, file_name);
        curr_im = rgb2gray(imread(curr_path));
        curr_im = histeq(curr_im, 255);
        curr_im = imresize(curr_im, im_size);
        curr_points = detectSURFFeatures(curr_im, ...
            'metricthreshold', 300, 'roi', im_roi);
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
                rays1, rays2, 0.999, thresh);
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
            [E_2pot, mask_2pot] = estimateRelativePose_PC2POT(...
                rays1, rays2, 0.999, thresh);
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
            disp([gt_rel(1:3, 1:3), normc(gt_rel(1:3, 4)), gt_rel(1:3, 4), gt_E]);
            disp('---');

%             return;
%             if i == 11
%                 return; 
%             end
            
        end
        [prev_points, prev_feat, prev_im] = deal(curr_points, curr_feat, curr_im);
        prev_i = i;
%         if sum(~isnan(t_err_5p)) > 50; return; end
    end
end


return;
