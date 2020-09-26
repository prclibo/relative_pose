
Q1 = rand(3,5);
Q2 = zeros(3,5);
V1 = zeros(3,5);
V2 = zeros(3,5);
lambda1 = zeros(1,5);
lambda2 = zeros(1,5);
U = zeros(3,5);

for k = 1:5
    Q1(:,k) = Q1(:,k)./norm(Q1(:,k),'fro');
    V1(:,k) = cross(rand(3,1),Q1(:,k));
    lambda1(k) = rand;
    U(:,k) = cross(V1(:,k),Q1(:,k))+lambda1(k)*Q1(:,k);
end

u0 = rand(4,1);
u0 = u0/norm(u0,'fro');
s = u0(4);
ss = s^2;
u = u0(1:3);
Rgt = (2*ss-1)*eye(3)+2*(u*u.'-s*skew(u));
Pgt = [Rgt rand(3,1)];

for k = 1:5
    p = Pgt*[U(:,k); 1];
    Q2(:,k) = p-U(:,k);
    Q2(:,k) = Q2(:,k)./norm(Q2(:,k),'fro');
    lambda2(k) = p.'*Q2(:,k);
    V2(:,k) = cross(Q2(:,k),U(:,k));
    %disp(Q1(:,k)'*V2(:,k)+Q2(:,k)'*V1(:,k));
end

Rt = g5pra(Q1,V1,Q2,V2,ss);

disp(reshape(Pgt,12,1));
disp(Rt);

