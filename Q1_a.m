clc;
clear all;
close all;

mu = [-2,2]; 
s1 = [1 0;0 1]; 
R1 = mvnrnd(mu,s1,1000);
a = -ones(1,height(R1)); 
scatter(R1(:,1),R1(:,2)) 
hold ON
mu = [2,-2];
s2 = [1 0;0 1];
R2 = mvnrnd(mu,s2,1000);
b = ones(1,height(R2));
w(1,:) = [1 1 -0.5]; 
scatter(R2(:,1),R2(:,2)) ;
x = -6:6;
y = (-w(3)-(w(1)*x))/w(2);

plot(x,y) %PLotting initial line
R = [R1;R2]; 
w0 = ones(height(R),1);
R = [R w0];
c = [a b];

%Obtaining initial efficiency
x1 = R(:,1);
x2 = R(:,2);

for i = 1:height(R)
    z1(i) = (w(1)*x1(i))+(w(2)*x2(i))+w(3);
    z(i) = z1(i)/abs(z1(i));
end
k = 1;
eff(k) = 1-(sum(z~=c)/length(z));
k1 = 0;%k1 will be changed from 0 in case of convergence


while k1 == 0
    k = k+1;
    for i = 1:height(R)
    z1(i) = (w(1)*x1(i))+(w(2)*x2(i))+w(3); 
    z(i) = z1(i)/abs(z1(i)); 
    end
    y = find(z~=c); 
    v = y(randperm(length(y),1));
    for i = 1: width(w)
        w(k,i) = w(k-1,i) +(0.5*(c(v)-z1(v))*R(v,i));%Modifying weight vectors
    end
    for i = 1:height(R)
    z1(i) = (w(k,1)*x1(i))+(w(k,2)*x2(i))+w(k,3);
    z(i) = z1(i)/abs(z1(i));
    end
    
    eff(k) =1 - (sum(z~=c)/length(z)); 
    if (eff(k) - eff(k-1) <= 0.1) 
        k1 = 1;
        if (eff(k) < eff(k-1))
            w(k,:)=w(k-1,:)
        end
    end
    

    if k1 == 0
        x = -6:6;
        y = (-w(k,3)-(w(k,1)*x))/w(k,2);
        plot(x,y)
    else
        x = -6:6;
        y = (-w(k,3)-(w(k,1)*x))/w(k,2);
        plot(x,y,'--r') %Plotting final line in dashed line
        title('Perceptron model with initial weight given');
    end
end
legend('class-1','class-2','boundary line-1','boundary line-2','final result')
hold OFF

