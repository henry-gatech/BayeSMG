%%% This script executes BayeSMG matrix completion on two image completion
%%% applications.
clearvars;
close all;

%% Lighthouse Image
load lighthouse;
[m1,m2] = size(X);
X = (X - mean(X(:)))/std(X(:));
X_org = X;
% configure how many entries we want to observe
prop = 0.5;
omega = zeros(m1,m2);
n = ceil(prop*m1*m2);
idx = randsample(m1*m2,n,false);
omega(idx) = 1; %set at sampled
% corrupt the matrix with noise
eta = 0.05;
for i=1:m1
    for j=1:m2
        X(i,j) = X(i,j) + randn*eta;
    end
end
% construct the partially observed matrix with noise
X = X.*omega;
% execute the BayeSMG completion method
r = 30;
[X_hat,lb,ub] = BayeSMG(X,omega,r,eta);
% plot the original matrix
figure
imagesc(X_org) 
colormap(hot(512))
caxis([-2.5 2.5])
colorbar
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
title('Original Matrix');
axis square
% plot the recovered matrix
figure
imagesc(X) 
colormap(hot(512))
caxis([-2.5 2.5])
colorbar
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
title('Partially Observed Matrix');
axis square

% plot the recovered matrix
figure
imagesc(X_hat) 
colormap(hot(512))
caxis([-2.5 2.5])
colorbar
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
title('Recovered Matrix by BayeSMG');
axis square
% plot the 95% HPD (UQ)
figure
imagesc(ub-lb)
colormap(jet)
colorbar
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
axis square
title('Uncertainty Interval (95% HPD) by BayeSMG');

%% Sunflare Image
img = imread('flare.tif');
X = double(img(1:256,1:256,1));
[m1,m2] = size(X);
X = (X - mean(X(:)))/std(X(:));
X_org = X;
% configure how many entries we want to observe
prop = 0.5;
omega = zeros(m1,m2);
n = ceil(prop*m1*m2);
idx = randsample(m1*m2,n,false);
omega(idx) = 1; %set at sampled
% corrupt the matrix with noise
eta = 0.05;
for i=1:m1
    for j=1:m2
        X(i,j) = X(i,j) + randn*eta;
    end
end
% construct the partially observed matrix with noise
X = X.*omega;
% execute the BayeSMG completion method
r = 30;
[X_hat,lb,ub] = BayeSMG(X,omega,r,eta);
% plot the original matrix
figure
imagesc(X_org) 
colormap(hot(512))
caxis([-2.5 2.5])
colorbar
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
title('Original Matrix');
axis square
% plot the recovered matrix
figure
imagesc(X) 
colormap(hot(512))
caxis([-2.5 2.5])
colorbar
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
title('Partially Observed Matrix');
axis square

% plot the recovered matrix
figure
imagesc(X_hat) 
colormap(hot(512))
caxis([-2.5 2.5])
colorbar
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
title('Recovered Matrix by BayeSMG');
axis square
% plot the 95% HPD (UQ)
figure
imagesc(ub-lb)
colormap(jet)
colorbar
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
axis square
title('Uncertainty Interval (95% HPD) by BayeSMG');
