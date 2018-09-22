%contrast enhancement
close all
clc, clear;
filter_int=@(n)1/(n^2)*ones(n);
L = 256;
I_original = imread('Testimage3.tif');
%mistake; this is not the average
%int_filter = 1/9*[1 1 1 1 1 ; 1 1 1 1 1 ; 1 1 1 1 1];
%% integration - blur
%int_filter_5x5 = 1/(5^2)*ones(5);
int_filter_3x3 = filter_int(3);
int_filter_5x5 = filter_int(5);
int_filter_7x7 = filter_int(7);
int_filter_9x9 = filter_int(9);

I_blur = conv2(I_original,int_filter_5x5,'same');
I_blur = uint8(I_blur);
%% correcting negatives
I_blur = correct_negatives(I_blur);
I_bin = I_blur;
% I_bin(I_bin<75)=0; %threshold to practially binarize the img
% I_bin(I_bin>=75)=255; %threshold to practially binarize the img
I_bin(I_bin < 125)=0; %threshold to practially binarize the img
I_bin(I_bin >= 125)=255; %threshold to practially binarize the img
%% remove spades / nums
I_rem = remove_spades_numbers(I_bin);
%% integration - blur
I_blur_rem = conv2(I_rem,int_filter_3x3,'same');
I_blur_rem = conv2(I_blur_rem ,int_filter_5x5,'same');
I_blur_rem = conv2(I_blur_rem ,int_filter_3x3,'same');
I_blur_rem = uint8(I_blur_rem);

%% 1st derivative
% x_1st_deriv = [0 1 0; 1 -4 1; 0 1 0];
% y_1st_deriv = [-1 0 1; -2 0 2; -1 0 1];
x_1st_deriv = [1 1; -1 -1];
y_1st_deriv = [1 -1; 1 -1];
I_1s_tder_x = conv2(I_blur_rem, x_1st_deriv ,'same');
I_1s_tder_y = conv2(I_blur_rem, y_1st_deriv ,'same');

x_1st_deriv = binarize_img(x_1st_deriv,200)
y_1st_deriv = binarize_img(y_1st_deriv,200)
%% 2nd derivative
sec_deriv = [0 1 0; 1 -4 1; 0 1 0];
I_sec_deriv = conv2(I_blur_rem, sec_deriv,'same');
I_sec_deriv = correct_negatives(I_sec_deriv) %corrects negative pixels
I_sec_deriv = binarize_img(I_sec_deriv,200)

figure
imshow([I_1s_tder_x, I_1s_tder_y, I_sec_deriv]);
[row, column] = find(I_sec_deriv > 100);
indeces = [row, column];
y_min = min(indeces(:,1));
y_max = max(indeces(:,1));
x_min = min(indeces(:,2));
x_max = max(indeces(:,2));

ROI_x = linspace(x_min-5,x_min+5,7);
ROI_y = linspace(y_min-5,y_max+5,7);

%% display figs
figure
imshow([I_original, I_blur, I_bin,I_rem, I_sec_deriv]);
title('1. original VS 2. blur 5x5 VS 3. binarized VS 4. spades/no removed from binarized VS 5. sec derivative with contrast enhancement');

rotate_img(I_rem,I_original);

function I_rem = remove_spades_numbers(I)
    [m,n] = size(I);
    counter_black = 0;
    counter_white = 0;
    change_from_white_to_black = 0;
    change_from_black_to_white = 0;
    l = m*n;
    
    white_index = find(I > 100 );
    for i=white_index(1):l
        if I(i-1) >= 100 && I(i) < 100
            change_from_white_to_black = 1;
            change_from_black_to_white = 0;
        elseif I(i-1) < 100 && I(i) >= 100
            change_from_white_to_black = 0;
            change_from_black_to_white  = 1;
        end
        % gets rid of black shapes within the card
        if I(i) < 100 && change_from_white_to_black
            counter_black = counter_black  + 1;
            %black_indeces = i;
        end
        
        if change_from_black_to_white
            if counter_black > 0 && counter_black < 100
                I(i-counter_black:i) = 255;
            end
            counter_black = 0;
        end
        
    end
    
    % gets rid of remaining white noise
    for i=2:l
        if I(i-1) >= 100 && I(i) < 100
            change_from_white_to_black = 1;
            change_from_black_to_white = 0;
        elseif I(i-1) < 100 && I(i) >= 100
            change_from_white_to_black = 0;
            change_from_black_to_white  = 1;
        end
    
        if I(i) > 100 && change_from_black_to_white
            counter_white = counter_white  + 1;
        end

        if change_from_white_to_black
            if counter_white > 0 && counter_white < 10
                I(i-counter_white:i) = 0;
            end
            counter_white = 0;
        end

    end
    I_rem = I;
end

function [img, unique_pixel_intensities, P_r] = contrast_enhancement(img)
    [m,n] = size(img);
    pixel_count_total = m*n;
    img = double(reshape(img, [1, pixel_count_total]));
    unique_pixel_intensities = double(unique(img));
    for i=1:length(unique_pixel_intensities)
       unique_indeces{i} = find(img == unique_pixel_intensities(i));
       P_r(i) = length(unique_indeces{i})/pixel_count_total;
       img(unique_indeces{i}) = 255*sum(P_r(1:i));
    end
    img = reshape(img, [m, n]);
end
function S = correct_negatives(r) %corrects negative pixels

    index_negatives = find(r < 0);
    if ~isempty(index_negatives)
        % index_negatives
        [m, n] = size(r);
        r = double(r);
        r_reshaped = reshape(r, [1,m*n]);
        min_r = min(r_reshaped);
        max_r = max(r_reshaped);
        S=uint8(255.*(r- min_r)./(max_r-min_r));
    else
        S = r;
    end
end

function S = neg_map(L, r)
    S=uint8((L-1)-r);
end
function I = binarize_img(I,threshold)
    I(I < threshold) = 0;
    I(I >= threshold) = 255;
end

function rotate_img(I_rem,I_original)
    figure
    imshow([I_rem, I_original]);
    %% Extract the Corner Indices
    [row,col] = find(I_rem == 255);
    indices = [row,col];
    
    mins = min(indices);
    maxs = max(indices);

    minY = mins(1);
    minX = mins(2);
    maxY = maxs(1);
    maxX = maxs(2);

    [cornerX, cornerY] = find(I_rem(minY,:) == 255);
    corner1 = [cornerY(1),minY];

    [cornerX, cornerY] = find(I_rem(maxY,:) == 255);
    corner2 = [cornerY(1),maxY];

    [cornerX, cornerY] = find(I_rem(:,minX) == 255);
    corner3 = [minX,cornerX(1)];

    [cornerX, cornerY] = find(I_rem(:,maxX) == 255);
    corner4 = [maxX,cornerX(1)];

    hold on
    plot(corner1(1),corner1(2),'rx', 'MarkerSize', 20);
    plot(corner2(1),corner2(2),'gx', 'MarkerSize', 20);
    plot(corner3(1),corner3(2),'go', 'MarkerSize', 20);
    plot(corner4(1),corner4(2),'ro', 'MarkerSize', 20);

    %% Rotate Image
    angle = atand(abs(corner4(2)-corner2(2))/(corner4(1)-corner2(1)));


    if(abs(corner4(1)-corner2(1)) < 220)
        imro = imrotate(I_original,-angle);
    else
        imro = imrotate(I_original,90-angle);
    end
    figure
    imshow(imro);
end