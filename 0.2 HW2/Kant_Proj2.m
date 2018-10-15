% By Mohak Kant and Alexandre Soares
% 2018 fall, Image Processing class at Texas Tech University
% Project 2

close all; clc, clear
%% original
filt_n =@(n).27*ones(n)./(n^2);
I = imread('Proj2.tif');
[m,n] = size(I);

I = double(I) - conv2(I,filt_n(3),'same');

figure
subplot(231);
imagesc(I);
title('original img');
axis image
%% fft
I_fft_shft = fftshift(fft2(I));
lg_abs_I_fft = abs(log(1+I_fft_shft));

subplot(232);
imagesc(lg_abs_I_fft);
title('Shifted and Scaled(log(1+abs)) FFT');
axis image
%% removing background 
mean_fft_lgbs = mean(lg_abs_I_fft (:));
std_fft_lgbs = std(lg_abs_I_fft (:));
threshold = mean_fft_lgbs + 4.55*std_fft_lgbs;
background_i_lo = find(lg_abs_I_fft   < threshold );
lg_abs_I_fft_thresh = (lg_abs_I_fft >= threshold);
I_fft_shft(background_i_lo) = 0;% I_fft_shft(background_i_lo)*.0001;

%lg_abs_I_fft = log(1+abs(I_fft_shft));
I_ifft1 = abs(ifft2(fftshift(I_fft_shft)));

subplot(233);
imagesc(I_ifft1);
title('ifft with non-uniform illumination');
axis image
%removing non=uniform illumination
%lg_abs_I_fft  = conv2(lg_abs_I_fft ,filt_n(1),'same'); 



c_val = I_fft_shft(205,273);
I_fft_shft(200:210,273) = 0;
I_fft_shft(205,267:280) = 0;
I_fft_shft(205,273) = c_val;

%saving it now without the pixels
lg_abs_I_fft = log(1+abs(I_fft_shft));

subplot(234);
imagesc(lg_abs_I_fft_thresh);
title('fft log-abs with threshold applied');
axis image

I_ifft = abs(ifft2(fftshift(I_fft_shft)));
subplot(235);
imagesc(I_ifft);
title('ifft final img');
axis image

subplot(236);
I_target = imread('Proj2_Output.tif');
imagesc(I_target); 
title('target img');
axis image
colormap gray