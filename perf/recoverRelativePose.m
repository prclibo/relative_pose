function poses = recoverRelativePose(Es, varargin)

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

poses = repmat(struct('R', [], 't', [], 'E', []), size(Es, 1) / 3 * 4, 1);
index = 1;
for i = 1:3:size(Es, 1)
    [U, ~, V] = svd(Es(i:i+2, :));
    sign_ = sign(det(U * V'));
    Rs = {};
    Rs{1} = U * W * V' * sign_;
    Rs{2} = U * W' * V' * sign_;
    t_x = U * Z * U' * sign_;
    t = [-t_x(2, 3); t_x(1, 3); -t_x(1, 2)];
%     disp([Rs{1}, Rs{2}]);
    for j = 1:2
        poses(index) = struct('R', Rs{j}, 't', t, 'E', Es(i:i+2, :));
        poses(index + 1) = struct('R', Rs{j}, 't', -t, 'E', Es(i:i+2, :));
        index = index + 2;
    end
end

best_pose = struct('trace_diff', -inf);
if ~isempty(parser.Results.NearestPose)
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
end

end
