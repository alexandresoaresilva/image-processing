close all; clc, clear
%% original
I = imread('Proj2.tif');
figure
subplot(241);
imagesc(I);
title('original img');
axis vis3d
%% fft
I_fft = fft2(I);
I_fft_shft = fftshift(I_fft);
lg_abs_I_fft = log(1+abs(I_fft_shft));
mean_fft_lgbs = mean(lg_abs_I_fft(:));
std_fft_lgbs = std(lg_abs_I_fft(:));

threshold = mean_fft_lgbs + 2*std_fft_lgbs;
background_i_lo = find(lg_abs_I_fft  <= threshold );
object_i_lo = find(lg_abs_I_fft  > threshold );
%lg_abs_I_fft(background_i_lo ) = 0;

m = 205;
n = 273;
c_va = lg_abs_I_fft(m,n);

subplot(242);
imagesc(lg_abs_I_fft);
title('fft log-abs');
axis vis3d 

%%manipulating fft2
I_fft_shft2 = I_fft_shft;
indeces_bkgrd = find( ~(I_fft_shft2 > (max(max(I_fft_shft2))-162)));
%I_fft_shft2(indeces_bkgrd) = I_fft_shft2(indeces_bkgrd)*.8;
% I_fft_shft2  = double(correct_negatives(exp_map(1,.25, I_fft_shft2)));
lg_abs_I_fft_2 = log(1+abs(I_fft_shft2));
lg_abs_I_fft_2 = lg_abs_I_fft_2+lg_abs_I_fft_2.^2+lg_abs_I_fft_2.^3;
background_i_lo2 = find(I_fft_shft2  <= threshold );
object_i_lo2 = find(I_fft_shft2  > threshold );
%lg_abs_I_fft(object_i_lo) = lg_abs_I_fft(object_i_lo);
%lg_abs_I_fft_2 = lg_abs_I_fft.^2
% lg_abs_I_fft((m-6):(m+6),n) = 0;
% lg_abs_I_fft(m,(n-6):(n+6)) = 0;
% lg_abs_I_fft(background_i_lo) = 0;
% lg_abs_I_fft(m,n) = c_va;
% SE = strel('diamond',1)
% I_er = imerode(lg_abs_I_fft ,SE);
subplot(243);
imagesc(lg_abs_I_fft_2);
title('fft log-abs - raised to power');
axis vis3d 
%% ifft
c_val = I_fft_shft(205,273);


% I_fft_shft(200:210,273) = 0;
% I_fft_shft(205,267:280) = 0;

% I_fft_shft((m-115):(m+11),n) = 0;
% I_fft_shft(m,(n-11):(n+11)) = 0;
% I_fft_shft(background_i_lo) = 0;
% I_fft_shft(205,273) = c_val;
%I_g = I_fft_shft; 


I_fft_shft(background_i_lo) = 0;
I_ifft = abs(ifft2(fftshift(I_fft_shft)));
%I_ifft_8bit = I_ifft;

subplot(244);
imagesc(I_ifft);
title('ifft final img');
axis vis3d


subplot(245);
I_target = imread('Proj2_Output.tif');
imagesc(I_target); 
title('target img');
axis vis3d


I_target_fft = fftshift(fft2(I_target));
I_target_fft_logabs = log(1 + abs(I_target_fft));
subplot(246);
imagesc(I_target_fft_logabs); 
title('target img fft - ');
axis vis3d
colormap gray
% 
% lg_abs_I_fft_center = lg_abs_I_fft(170:250,250:300);
% lg_abs_I_fft_center2 = conv2(lg_abs_I_fft_center,norm_filt(7),'same');
% lg_abs_I_fft(171:249,251:299) = lg_abs_I_fft_center2(2:end-1,2:end-1);
% 
% I_ifft_lg_abs_I_fft = abs(ifft2(fftshift(lg_abs_I_fft.^1.5)));

% figure
% imagesc(lg_abs_I_fft_center); 
% title('Center');
% axis vis3d
% colormap gray
% 
% figure
% imagesc(lg_abs_I_fft_center2); 
% title('Blurred Center');
% axis vis3d
% colormap gray
% 
% figure
% imagesc(lg_abs_I_fft); 
% title('Blurred Center Replaced');
% axis vis3d
% colormap gray

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
            filt(i,j) = guass_filt(x(i),y(j),1);
        end
    end
end

function S = exp_map(C, alpha, r) %exponential mapping
    %S=C.*exp(double(r));
    r = double(r);
    S=C.*r.^(alpha);
end