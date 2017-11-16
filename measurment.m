clc
close all

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

for k = 1:size_w
    theta(k) = rand(1)*2*pi;
end 

alpha = (1/sqrt(2))*[1; 1];
%theta = 0:0.01:2*pi;
w = ro*exp(i*theta);

for k = 1:size_w
    a = (s*alpha + (1/sqrt(2))*w(k)*[1; -1]);
    %a = (1/sqrt(2))*[s + w(k); s - w(k)];
    a_mod1(k) = abs(a(1));
    a_mod2(k) = abs(a(2));
    gamma(k) = 1;
end 

scatter(theta, abs(a_mod1 - a_mod2))
figure

scatter(theta, a_mod1)
hold on
scatter(theta, a_mod2)

% plot(theta, gamma);
% hold on
% plot(theta, a_mod1, 'r', theta, a_mod2, 'b');
% xlabel('theta (rad)');
% ylabel(' ')

%% Another Two-Dimensional Example
%% Para Qubit |0>
clc
close all

size_w = 158;
ro = 1;
gamma = ro;
s = (sqrt(2)-1)*ro;
% for k = 1:size_w
%     theta(k) = rand(1)*(pi/2);
% end

% for k = 1:size_w
%     psi_1(k) = rand(1)*(2*pi);
% end
% 
% for k = 1:size_w
%     psi_2(k) = rand(1)*(2*pi);
% end

psi_1 = 0;
psi_2 = 0;

theta = 0:0.01:pi/2;
%psi_1 = 0:0.01:2*pi;
%psi_2 = 0:0.01:2*pi;

alpha = [1;0];
w = ro*[exp(i*psi_1)*cos(theta); exp(i*psi_2)*sin(theta)];

for k = 1:size_w
    w(k)
    a = s*alpha + ro*[exp(i*psi_1)*cos(theta(k)); exp(i*psi_2)*sin(theta(k))];
    a_mod1(k) = abs(a(1));
    a_mod2(k) = abs(a(2));
    gamma(k) = ro;
end 

% scatter(theta, a_mod1 - a_mod2)
% figure
% scatter(theta, a_mod1)
% hold on
% scatter(theta, a_mod2)
plot(theta, gamma);
hold on
plot(theta, a_mod1, 'r', theta, a_mod2, 'b');
xlabel('theta (rad)');
ylabel(' ')

%% para Qubit = 1/sqrt(2) * ( |0> + |1> )

clc
close all

size_w = 158;
ro = 1;
gamma = ro;
s = (sqrt(2)-1)*ro;
% for k = 1:size_w
%     theta(k) = rand(1)*(pi/2);
% end

% for k = 1:size_w
%     psi_1(k) = rand(1)*(2*pi);
% end
% 
% for k = 1:size_w
%     psi_2(k) = rand(1)*(2*pi);
% end

psi_1 = 0;
psi_2 = 0;

theta = 0:0.01:pi/2;
%psi_1 = 0:0.01:2*pi;
%psi_2 = 0:0.01:2*pi;

alpha = [1;0];
w = ro*[exp(i*psi_1)*cos(theta); exp(i*psi_2)*sin(theta)];

for k = 1:size_w
    w(k)
    a = s*alpha + ro*[exp(i*psi_1)*cos(theta(k)); exp(i*psi_2)*sin(theta(k))];
    a_mod1(k) = abs(a(1))^2;
    a_mod2(k) = abs(a(2))^2;
    gamma(k) = ro;
end 

% scatter(theta, a_mod1 - a_mod2)
% figure
% scatter(theta, a_mod1)
% hold on
% scatter(theta, a_mod2)
plot(theta, gamma);
hold on
plot(theta, a_mod1, 'r', theta, a_mod2, 'b');
xlabel('theta (rad)');
ylabel(' ')

%% para Qubit = a*|0> + b*|1>

clc
close all

size_w = 158;
ro = 1;
gamma = ro;
s = (sqrt(2)-1)*ro;
% for k = 1:size_w
%     theta(k) = rand(1)*(pi/2);
% end

% for k = 1:size_w
%     psi_1(k) = rand(1)*(2*pi);
% end
% 
% for k = 1:size_w
%     psi_2(k) = rand(1)*(2*pi);
% end

psi_1 = 0;
psi_2 = 0;

theta = 0:0.01:pi/2;
%psi_1 = 0:0.01:2*pi;
%psi_2 = 0:0.01:2*pi;

q = 0:0.01:1;

w = ro*[exp(i*psi_1)*cos(theta); exp(i*psi_2)*sin(theta)];

plot(theta, gamma);
hold on

for l = 1:101
    alpha = [q(l);1-q(l)];
    for k = 1:size_w
        a = s*alpha + ro*[exp(i*psi_1)*cos(theta(k)); exp(i*psi_2)*sin(theta(k))];
        a_mod1(k) = abs(a(1))^2;
        a_mod2(k) = abs(a(2))^2;
        gamma(k) = ro;
    end 
    plot(theta, a_mod1, 'r', theta, a_mod2, 'b');
    xlabel('theta (rad)');
    ylabel(' ')
    hold on
end

% scatter(theta, a_mod1 - a_mod2)
% figure
% scatter(theta, a_mod1)
% hold on
% scatter(theta, a_mod2)

