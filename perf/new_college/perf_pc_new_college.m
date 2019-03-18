data_dir = '/Users/li/data/new_college/CityCentre';
im_dir = fullfile(data_dir, 'Images');
pos_path = fullfile(data_dir, 'ImageCollectionCoordinates.mat');

addpath(fullfile(fileparts(mfilename('fullpath')), '/../../build/matlab/'));
addpath(fullfile(fileparts(mfilename('fullpath')), '/../utils'));

load(pos_path);
t = cumsum([0; sqrt(sum(diff(GPS) .^ 2, 2))]);

splx = spline(t, GPS(:, 1));
sply = spline(t, GPS(:, 2));
tang = normr([ppval(fnder(sply),t), ppval(fnder(splx),t)]);

im_indices = intersect(find(isfinite(tang)), 50:2:size(GPS, 1));

virt_extrs = [1, 0, 0, 0; 0, 0, 1, 0; 0, -1, 0, 0; 0, 0, 0, 1];

min_move = 1;
for index = 1:numel(im_indices)
    i = im_indices(index);
    fprintf('%d / %d\n', i, numel(im_indices));
    rng(3);
%     clc
    if index > 1
        curr_pose = [tang(i, 2), tang(i, 1), 0, GPS(i, 1) - GPS(prev_i, 1);
                    -tang(i, 1), tang(i, 2), 0, GPS(i, 2) - GPS(prev_i, 2);
                    0, 0, 1, 0; 0, 0, 0, 1];
        prev_pose = [tang(prev_i, 2), tang(prev_i, 1), 0, 0;
                    -tang(prev_i, 1), tang(prev_i, 2), 0, 0;
                    0, 0, 1, 0; 0, 0, 0, 1];
        gps_rel = prev_pose \ curr_pose;
        gps_move = norm(gps_rel(1:3, 4));
%         gt_rel = virt_extrs \ gps_rel * virt_extrs;
%         
%         gt_nt = normc(gt_rel(1:3, 4));
%         gt_move = norm(gt_rel(1:3, 4));
%         gt_E = skew(gt_nt) * gt_rel(1:3, 1:3);
%         gt_E = gt_E / norm(gt_E(:));
%         disp([gt_rel(1:3, 1:3), normc(gt_rel(1:3, 4)), gt_rel(1:3, 4), gt_E]);            
    end
    if index == 1 || gps_move > min_move
        file_name = sprintf('%04d.jpg', i);
        curr_path = fullfile(im_dir, file_name);
        curr_im = rgb2gray(imread(curr_path));
        curr_im = histeq(curr_im, 255);
        curr_points = detectSURFFeatures(curr_im, ...
            'metricthreshold', 300); % , 'roi', im_roi);
        [curr_feat, curr_points] = extractFeatures(curr_im, curr_points);
        if index > 1% && gt_move > min_move
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
            end
        end
    end
    [prev_points, prev_feat, prev_im] = deal(curr_points, curr_feat, curr_im);
    prev_i = i;
end