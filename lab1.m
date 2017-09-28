%% lab 1

% 2

clc
clf
clear all 
close all

N = 2000;
x_k = randn([N,1]);          % Normal
u_k = rand([N,1]);           % Uniform

% hist(x_k,10)
% title('Standard Normal Distribution')
% xlabel('x_k')
% ylabel('Frequency')
% figure
% hist(u_k,10)
% title('Standard Uniform Distribution')
% xlabel('u_k')
% ylabel('Frequency')
% figure

X = randn([N,1]);               % Normal
Y = randn([N,1]);               % Normal 
%Y = exprnd(0,[N,1]);


a = sqrt(3);
u = -a + 2*a*rand([N,1]);       % Uniform (-a,a)
v = -a + 2*a*rand([N,1]);       % Uniform (-a,a)

% plot(X,Y,'.')
% title('Gaussian vs. Gaussian')    % Gaussian is centered around its mean which in this case is 0. This can be seen from the histogram. It's more likely to find a point in the region close to zero. 
% xlabel('x_k')
% ylabel('y_k')
% figure
% plot(u,v,'.')
% title('Uniform vs. Uniform')      % The points are distributed uniformly which can be seen from the histogram. It's equally likely to find a point in the center as in one of the edges. 
% xlabel('u_k')
% ylabel('v_k')
% figure
% plot(X,u,'.')
% title('Gaussian vs. Uniform')     % The Gaussian r.v is on the x-axis which means the points will be distributed around x=0 (the mean of the Gaussian). 
% xlabel('x_k')
% ylabel('u_k')




Y_hat = 0.5;     % Fixed value of Y
dY = 1.1;        % Tolerance 

minY = Y_hat - dY;
maxY = Y_hat + dY;


Xk = zeros(N,1);
for i=1:N
if Y(i) > minY && Y(i) < maxY
Xk(i,1) = X(i);
end

end
Xk(Xk==0) = [];
histogram(Xk,10)
title('Conditional Distribution')
xlabel('$$x|y=\hat{y}$$','Interpreter','Latex')
ylabel('Frequency')
figure
histogram(X,10)
title('Marginal distribution')
xlabel('x')
ylabel('Frequency') 


%% 3.1

clc
clf
clear all 
close all

N = 2000;
alpha = [0.5 -0.5 0.9 -0.9];
X = randn([N,1]);
Y = randn([N,1]);




Z1 = alpha(1) * X + sqrt(1-alpha(1).^2) * Y; 
Z2 = alpha(2) * X + sqrt(1-alpha(2).^2) * Y; 
Z3 = alpha(3) * X + sqrt(1-alpha(3).^2) * Y; 
Z4 = alpha(4) * X + sqrt(1-alpha(4).^2) * Y; 

% subplot(2,2,1);
% plot(X,Z1,'.')
% title('\alpha=0.5')
% xlabel('x')
% ylabel('z')
% hold on
% plot(X,alpha(1)*X)            % r_xz = alpha and r_yz = sqrt(1-alpha)^2
% 
% subplot(2,2,2);
% plot(X,Z2,'.')
% title('\alpha=-0.5')
% xlabel('x')
% ylabel('z')
% hold on
% plot(X,alpha(2)*X)            % r_xz = alpha and r_yz = sqrt(1-alpha)^2
% 
% subplot(2,2,3);
% plot(X,Z3,'.')
% title('\alpha=0.9')
% xlabel('x')
% ylabel('z')
% hold on
% plot(X,alpha(3)*X)            % r_xz = alpha and r_yz = sqrt(1-alpha)^2
% 
% 
% subplot(2,2,4);
% plot(X,Z4,'.')
% title('\alpha=-0.9')
% xlabel('x')
% ylabel('z')
% hold on
% plot(X,alpha(4)*X)            % r_xz = alpha and r_yz = sqrt(1-alpha)^2


% figure;
% histogram(Z3)
% title('\alpha = 0.9')
% xlabel('z')
% ylabel('Frequency')


%% 3.2
clc
clf
clear all 
close all

N = 2000;
alpha = [0.7 -0.7];
X = randn([N,1]);
Y = randn([N,1]);




Z1 = alpha(1) * X + sqrt(1-alpha(1).^2) * Y; 
Z2 = alpha(2) * X + sqrt(1-alpha(2).^2) * Y; 


Z_hat1 = -0.5;     % Fixed value of Z
Z_hat2 = 0.5;
dZ = 0.1;        % Tolerance 

minZ1 = Z_hat1 - dZ;
maxZ1 = Z_hat1 + dZ;

minZ2 = Z_hat2 - dZ;
maxZ2 = Z_hat2 + dZ;

Xk1 = zeros(N,1);
for i=1:N
    if Z1(i) > minZ1 && Z1(i) < maxZ1
        Xk1(i,1) = X(i);
    end
    
end

Xk2 = zeros(N,1);
for i=1:N
    if Z1(i) > minZ2 && Z1(i) < maxZ2
        Xk2(i,1) = X(i);
    end
    
end

Xk3 = zeros(N,1);
for i=1:N
    if Z2(i) > minZ1 && Z2(i) < maxZ1
        Xk3(i,1) = X(i);
    end
    
end

Xk4 = zeros(N,1);
for i=1:N
    if Z2(i) > minZ2 && Z2(i) < maxZ2
        Xk4(i,1) = X(i);
    end
    
end
Xk1(Xk1==0) = [];
Xk2(Xk2==0) = [];
Xk3(Xk3==0) = [];
Xk4(Xk4==0) = [];
binlimits = [-2 2];
subplot(2,2,1)
histogram(Xk1,'BinLimits',binlimits)               % The mean of Xk is given by the condition on Z. 
title('$$\alpha = 0.7$$','Interpreter','Latex')
xlabel('$$x|z=-0.5$$','Interpreter','Latex')
ylabel('Frequency')
subplot(2,2,2)
histogram(Xk2,'BinLimits',binlimits)               % The mean of Xk is given by the condition on Z. 
title('$$\alpha = 0.7$$','Interpreter','Latex')
xlabel('$$x|z=0.5$$','Interpreter','Latex')
ylabel('Frequency')
subplot(2,2,3)
histogram(Xk3,'BinLimits',binlimits)               % The mean of Xk is given by the condition on Z. 
title('$$\alpha = -0.7$$','Interpreter','Latex')
xlabel('$$x|z=-0.5$$','Interpreter','Latex')
ylabel('Frequency')
subplot(2,2,4)
histogram(Xk4,'BinLimits',binlimits)               % The mean of Xk is given by the condition on Z. 
title('$$\alpha = -0.7$$','Interpreter','Latex')
xlabel('$$x|z=0.5$$','Interpreter','Latex')
ylabel('Frequency')



%% 4 
clc
clf
clear all
close all

N = 256;
K = N;
X = randn(N,K);

eAverage = mean(X,2);     % Mean of all rows
tAverage = mean(X,1);     % Mean of all columns

% plot(eAverage)
% title('Time Average vs. Ensemble Average')
% hold on
% plot(tAverage)
% axis([0 256 -0.2 0.2])
% xlabel('Samples and time realization')
% ylabel('Mean')

% ERGODIC!!!!
n1 = 3;
n2 = 2;

x1 = X(n1,:);
x2 = X(n2,:);
plot(x1,x2,'.')
title('Joint Distribution')
xlabel('x[n_1]')
ylabel('x[n_2]')

sum = 0;
for i=1:K
    product = X(n1,i).*X(n2,i);
    sum = sum + product;   
end

sum = sum / K;

% UNCORRELATED!!!!!!!!!!!!!!


%% 5 
clc
clf
clear all
close all



