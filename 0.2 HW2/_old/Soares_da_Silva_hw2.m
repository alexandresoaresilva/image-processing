close all; clc, clear

filt_n = @(n)1/n^2*ones(n);
adjust_lg_abs_I_fft =...
    @(lg_abs_I_fft)uint8(round(lg_abs_I_fft./max(max(lg_abs_I_fft))*255));
I_read = imread('Proj2.tif');
I = I_read;
%I_conv = conv2(I,filt_n(9),'same');
%I_conv = I_read - uint8(I_conv);
%I = correct_negatives(I_conv)

I_fft = fft2(I);
abs_I_fft = abs(fftshift(I_fft));
lg_abs_I_fft = log(abs_I_fft);
figure


pattern_i = find(lg_abs_I_fft > 12 );
background_i_lo = find(lg_abs_I_fft <= 14);
I_g = fftshift(I_fft); 
I_g(background_i_lo ) = 0;
% I_g(background_i_hi) = 0;

%I_ifft = ifft2(fftshift(I_g));
I_ifft = ifft2(fftshift(I_g));
I_ifft = abs(I_ifft);
% a = uint8(round(I_ifft/max(max(I_ifft))*255));
a = uint8(round(I_ifft));
[m,n] = size(lg_abs_I_fft);

a = histeq(a);
imagesc([I, lg_abs_I_fft, lg_abs_I_fft; lg_abs_I_fft, I_g, a]);

%imagesc([I, lg_abs_I_fft, lg_abs_I_fft; lg_abs_I_fft, I_g, a]);
title('original, lg_abs_I_fft, I_conv');
colormap gray

function S = correct_negatives(r) %corrects negative pixels
    [m, n] = size(r);
    r = double(r);
    r_reshaped = reshape(r, [1,m*n]);
    
    min_r = min(r_reshaped);
    max_r = max(r_reshaped);
    S=uint8(255.*(r- min_r)./(max_r-min_r));
end

function binImage =  binI(threshold, image)
    binImage = image;
    binImage(image<(threshold)) = 0;
    binImage(image>(threshold)) = 255;
end

function filt = norm_filt(n)
    guass_filt =@(x,y,std)exp(-(x.^2+y^2)/(2*std^2))*1/(2*pi*std^2);

    x = -n:1:n;
    y= x;
    len = length(x);
    for i=1:len
        for j=1:len
            filt(i,j) = guass_filt(x(i),y(j),4);
        end
    end
end