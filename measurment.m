clc
close all

1

%%  Simple Two-Dimensional Example
%% para alpha = |0> ; w = ro*exp(i*theta)*[0; 1]

close all

ro = 1;
s = 1*ro;  

theta = 0:0.01:2*pi;

alpha = s*[1;0];
w = ro*exp(i.*theta);

for k = 1:size_w
     a = s*alpha +  w(k)*[0; 1];
     a_mod1(k) = abs(a(1))^2;
     a_mod2(k) = abs(a(2))^2;
end 

plot(theta, a_mod1, 'r', theta, a_mod2, 'b')


%% Para alpha = |0> + |1> ; s > ro ; ro = gamma ; w = ro*exp(i*theta)*[1; -1]

close all

size_w = 629;
ro = 1;
s = 1;

%for k = 1:size_w
%    theta(k) = rand(1)*pi;
%end 

alpha = (1/sqrt(2))*[1; 1];
theta = 0:0.01:2*pi;
w = ro*exp(i*theta);

for k = 1:size_w
    a = (s*alpha + (1/sqrt(2))*w(k)*[1; -1]);
    %a = (1/sqrt(2))*[s + w(k); s - w(k)];
    a_mod1(k) = abs(a(1));
    a_mod2(k) = abs(a(2));
    gamma(k) = 1;
end 
plot(theta, a_mod1+a_mod2)
figure

plot(theta, gamma);
hold on
plot(theta, a_mod1, 'r', theta, a_mod2, 'b');
xlabel('theta (rad)');
ylabel(' ')
