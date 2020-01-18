function d1wlsdestripedemo
%D1WLSDESTRIPEDEMO CNLM denoising method.
%  D1WLSDESTRIPEDEMO reads and denoises raw uncooled infrared images with
%  stripe noise.
%
%  To run the demo, type D1WLSDESTRIPEDEMO from the Matlab prompt.

%  Fangzhou Li
%
%  November 2019

disp(' ');
disp('  **********  D1-WLS Destriping Demo  **********');
disp(' ');
disp('  This demo reads raw uncooled infrared images with stripe noise.');
disp('  The stripe noise will be removed using two 1D filters on horizontal and vertical direction.');
disp('  The denoised image and the estimated stripe noise will be shown.');
disp(' ');

%% prompt user for image %%

addpath('.\functions')

im = readImage('d1wlsdestripedemo');

%% generate noisy image %%

% sigma_white = 0;
% sigma_stripe = 0;

% disp(' ');
% disp('Generating mixed noisy image...');

% n = randn(size(im)) * sigma_white;
% n = n + repmat(randn(1, size(im, 2)), size(im, 1), 1) .* sigma_stripe;
% imnoise = im + n;

if length(size(im)) > 2
    im = double(rgb2gray(im));
else
    im = double(im);
end    

% denoise!
disp('Performing D1WLS destriping...');
[dI] = d1_WLS_Destriping(im, 40, 3);

% show results %

figure; imshow(newlp(im));
title('Original noisy image');

figure; imshow(newlp(dI));
title('Denoised image')

figure; imshow(newlp(im - dI))
title('Extracted stripe noise')
