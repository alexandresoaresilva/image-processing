close all
clc, clear;
r = rgb2gray(imread('paris_879.jpg'));

subplot(141)
imshow(r);
axis vis3d 

L = 256;
subplot(142)
imshow(neg_map(L, r));
axis vis3d 
C= 1;
subplot(143)
imshow(correct_negatives(log_map(C, r)));
axis vis3d 
C= 1;
alpha = 0.8;
subplot(144)
imshow(correct_negatives(exp_map(C,alpha, r)));
axis vis3d 
%close all

%negative mapping
function S = neg_map(L, r)
    S=uint8((L-1)-r);
end

function S = log_map(C, r) %log mapping
    S=C*log(1 + double(r));
end

function S = exp_map(C, alpha, r) %exponential mapping
    %S=C.*exp(double(r));
    r = double(r);
    S=C.*r.^(alpha);
end

function S = correct_negatives(r) %corrects negative pixels
    [m, n] = size(r);
    r = double(r);
    r_reshaped = reshape(r, [1,m*n]);
    
    min_r = min(r_reshaped);
    max_r = max(r_reshaped);
    S=uint8(255.*(r- min_r)./(max_r-min_r));
end