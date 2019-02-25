clear;

sdk_dir = '~/workspace/robotcar-dataset-sdk/matlab';
image_dir = '~/data/robotcar/sample/stereo/centre';
timestamp_dir = '~/data/robotcar/sample/stereo.timestamps';
model_dir = '~/workspace/robotcar-dataset-sdk/models';
ins_file = '~/data/robotcar/sample/gps/ins.csv';
consec_interv = 5;

addpath(sdk_dir);
timestamps = load(timestamp_dir);
timestamps = timestamps(1:consec_interv:end, 1)';
[fx, fy, cx, cy, G_camera_image, LUT] = ReadCameraModel(image_dir, model_dir);
otime = timestamps(1);

poses = InterpolatePoses(ins_file, timestamps, otime);

for i = 1:numel(timestamps)
    curr_im = rgb2gray(LoadImage(image_dir, timestamps(i), LUT));
    curr_points = detectSURFFeatures(curr_im);
    [curr_feat, curr_points] = extractFeatures(curr_im, curr_points);
    
    if i > 1
        pairs = matchFeatures(curr_feat, prev_feat, 'method', 'approximate') ;
        mpoints1 = curr_points(pairs(:, 1));
        mpoints2 = prev_points(pairs(:, 2));
        rays1 = ([mpoints1(:).Location] - [cx, cy]) ./ [fx, fy];
        rays2 = ([mpoints2(:).Location] - [cx, cy]) ./ [fx, fy];
        [rays1(:, 3), rays2(:, 3)] = deal(1, 1);
        
        
        figure; showMatchedFeatures(curr_im, prev_im, mpoints1, mpoints2);

    end
    
    [prev_points, prev_feat, prev_im] = deal(curr_points, curr_feat, curr_im);
end