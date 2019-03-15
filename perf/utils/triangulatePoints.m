function [trgl1, trgl2, mag] = triangulatePoints(rays1, rays2, pose)

trgl1 = nan(size(rays1));
mag = nan(size(rays1, 1), 1);
for i = 1:size(rays1, 1)
    A1 = skew(rays1(i, :)) * eye(3, 4);
    A2 = skew(rays2(i, :)) * [pose.R, pose.t];
    [~, S, V] = svd([A1; A2]);
%     X = [A1; A2] \ zeros(6, 1);
    trgl1(i, :) = V(1:3, 4) ./ V(4, 4);
    mag(i) = S(4, 4);
end

trgl2 = bsxfun(@plus, trgl1 * pose.R', pose.t');

end
