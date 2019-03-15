function calib = loadCalib(calib_path)

fid = fopen(calib_path, 'r');
calib_list = textscan(fid, '%s %f %f %f %f %f %f %f %f %f %f %f %f', ...
    'Delimiter', {' ', ': '});

fclose(fid);

calib = repmat(struct('extrs', [], 'intrs', []), 4, 1);
for i = 1:4
    P = nan(12, 1);
    for j = 1:12; P(j) = calib_list{j + 1}(i); end
    P = reshape(P, 4, 3)';
    [R, K] = qr(flip(P(1:3, 1:3))');
    K = rot90(K', 2);
    R = fliplr(R');
    S = diag(sign(diag(K)));
    K = K * S;
    R = S * R;
    t = R' * P(:, 4);
    calib(i).extrs = [R, t; 0, 0, 0, 1];
    calib(i).intrs = K;
end