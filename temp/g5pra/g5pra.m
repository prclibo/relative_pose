% 5-point solver, generalized camera and known relative angle
function Rt = g5pra(Q1,V1,Q2,V2,ss)

    s = sqrt(ss);
    t = sqrt(1-ss);

    B = matrix37x81(Q1,V1,Q2,V2,ss,s);

    B = B(:,1:37)\B(:,38:81);

    % 44 x 44 action matrix
    C = [-B([17,26,27,28,29,30,31], :); zeros(1,44); -B([34,35,36,37],:); zeros(32,44)];
    C(8,1) = 1;
    C(13,2) = 1;
    C(14,3) = 1;
    C(15,4) = 1;
    C(16,5) = 1;
    C(17,6) = 1;
    C(18,7) = 1;
    C(19,8) = 1;
    C(20,11) = 1;
    C(21,12) = 1;
    C(22,13) = 1;
    C(23,14) = 1;
    C(24,15) = 1;
    C(25,16) = 1;
    C(26,17) = 1;
    C(27,18) = 1;
    C(28,19) = 1;
    C(29,22) = 1;
    C(30,23) = 1;
    C(31,24) = 1;
    C(32,25) = 1;
    C(33,26) = 1;
    C(34,27) = 1;
    C(35,28) = 1;
    C(36,31) = 1;
    C(37,32) = 1;
    C(38,33) = 1;
    C(39,34) = 1;
    C(40,35) = 1;
    C(41,38) = 1;
    C(42,39) = 1;
    C(43,40) = 1;
    C(44,43) = 1;

    % eigenvectors of action matrix encode solutions
    [V,~] = eig(C);

    Axs = V(41:44,:);
    Axs = Axs(:, imag(Axs(1,:))==0 & imag(Axs(2,:))==0 & imag(Axs(3,:))==0 & Axs(4,:)~=0);
    nsr = size(Axs,2); % number of real solutions
    Rt = zeros(12,nsr);

    for j = 1:nsr
        u = Axs(1:3,j)./Axs(4,j);
        u = (u./sqrt(u.'*u)).*t;
        Rt(:,j) = Rt_g5pra(Q1,V1,Q2,V2,u,ss,s);
    end

end