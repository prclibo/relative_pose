
function rt = Rt_g5pra(Q1,V1,Q2,V2,u,ss,s)
    % rotation matrix R
    R = (2*ss-1)*eye(3)+2*(u*u.'-s*skew(u));

    S = zeros(5,4);
    for i = 1:5
        S(i,1:3) = -cross(R*Q1(:,i),Q2(:,i));
        S(i,4) = Q2(:,i).'*R*V1(:,i)+Q1(:,i).'*R.'*V2(:,i);
    end

    % translation vector t
    [~,~,V] = svd(S,0);
    t0 = V(:,4);
    t = t0(1:3)./t0(4);
    rt = [reshape(R,9,1); t];
end