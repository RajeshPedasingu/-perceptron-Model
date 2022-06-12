clc;
clear all;
close all;

mu1 = [-2,2];
s = [1 0;0 1];
R1 = mvnrnd(mu1,s,1000);
a = -ones(1,height(R1));

mu2 = [2,-2];
s = [1 0;0 1];
R2 = mvnrnd(mu2,s,1000);
b = ones(1,height(R2));
figure()
scatter(R1(:,1),R1(:,2))
hold ON
scatter(R2(:,1),R2(:,2))
title('Scatter diagram before Baysian classification directly from random generation');
hold OFF

c = [a b];
R = [R1;R2];
d = [-1,1];
for i = 1:height(R)
    S = [mvnpdf(R(i,:),mu1,s),mvnpdf(R(i,:),mu2,s)];
    S1(i) = find(S==max(S));
    M(i) = d(S1(i));
end
eff = 1 - (sum(M~=c)/length(M))
R1 = R(M==-1,:);
R2 = R(M==1,:);
figure()
scatter(R1(:,1),R1(:,2))
hold ON
scatter(R2(:,1),R2(:,2))
title('Scatter diagram after Baysian classification');
hold OFF

