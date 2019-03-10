function poses = recoverRelativePose(Es, varargin)

MAX_DEPTH = 50;

parser = inputParser();
parser.StructExpand = false;
parser.addRequired('Es', @ismatrix);
parser.addOptional('rays1', [], @isnumeric);
parser.addOptional('rays2', [], @isnumeric);
parser.addParameter('angle', [], @isscalar);
parser.addParameter('ZeroScrewTransl', false, @islogical);
parser.addParameter('threshold', 5e-3, @isscalar);
parser.addParameter('NearestPose', struct([]), @isstruct);
parser.parse(Es, varargin{:});

W = [0, -1, 0; 1, 0, 0; 0, 0, 1];
Z = [0, -1, 0; 1, 0, 0; 0, 0, 0];

poses = repmat(struct('R', [], 't', [], 'E', []), 0, 1);
index = 1;
for i = 1:3:size(Es, 1)
    [U, ~, V] = svd(Es(i:i+2, :));
    sign_ = sign(det(U * V'));
    Rs = {};
    Rs{1} = U * W * V' * sign_;
    Rs{2} = U * W' * V' * sign_;
    t_x = U * Z * U' * sign_;
    t = [-t_x(2, 3); t_x(1, 3); -t_x(1, 2)];
    t = normc(t);
%     disp([Rs{1}, Rs{2}]);
    for j = 1:2
        if parser.Results.ZeroScrewTransl
            rotvec = vrrotmat2vec(Rs{j});
            if abs(rotvec(end) - pi) < 1e-4; continue; end
        end
        poses(index) = struct('R', Rs{j}, 't', t, 'E', Es(i:i+2, :));
        poses(index + 1) = struct('R', Rs{j}, 't', -t, 'E', Es(i:i+2, :));
        index = index + 2;
    end
end

if ~isempty(parser.Results.NearestPose)
    best_pose = struct('trace_diff', -inf);
    % Omit every other mirrored pose
    for i = 1:2:numel(poses)
        relative = poses(i).R' * parser.Results.NearestPose.R;
        if trace(relative) > best_pose.trace_diff
            best_pose = poses(i);
            best_pose.i = (i + 1) / 2;
            best_pose.trace_diff = trace(relative);
        end
    end
    poses = best_pose;
    cosine = max(-1, min(1, (poses.trace_diff - 1) / 2));
    poses.rot_diff = acosd(cosine);
    cosine = max(-1, min(1, dot(poses.t, parser.Results.NearestPose.t)));
    transl_diff = acosd(cosine);
    poses.transl_diff = min(transl_diff, 180 - transl_diff);
    return;
elseif ~isempty(parser.Results.rays1)
    best_pose = struct('pos_depth', 0);
    for i = 1:numel(poses)
        [trgl1, trgl2] = triangulatePoints(...
            parser.Results.rays1, parser.Results.rays2, poses(i));
        pos_depth = sum(trgl1(:, 3) > 0 & trgl1(:, 3) < MAX_DEPTH &...
            trgl2(:, 3) > 0 & trgl2(:, 3) < MAX_DEPTH);
        if pos_depth > best_pose.pos_depth
            best_pose = poses(i);
            best_pose.pos_depth = pos_depth;
        end
    end
    poses = best_pose;
end

end

