%written by: Alexandre Soares & Mohak Kant
% HW for ECE 4367 - Image Processing at Texas Tech University, prof. Dr Sari-Sarraf (fall 2018)
% Sep 27, 2018

close all
clc, clear;
filter_int=@(n)1/(n^2)*ones(n);
L = 256;
I_original = imread('Testimage5.tif');
%% integration - blur
I_blur = conv2(I_original,filter_int(13),'same');
I_blur = uint8(I_blur);
%% correcting negatives
I_blur = correct_negatives(I_blur); %just if there are negatives
I_bin = binarize_img(I_blur,125);
%% remove spades / nums and rotate
I_bin_no_spades = remove_spades_numbers(I_bin);
[I_rotat, I_rotated_bin] = rotate_img(I_bin_no_spades,I_original);
I_cropped = crop_img(I_rotated_bin, I_rotat);
I = check_if_upside_down(I_cropped);
%% display figs
figure
hold on
subplot(221);
imshow([I_original, I_blur]);
title('Original vs 5x5 blur (int)');
subplot(222);
imshow([I_bin, I_bin_no_spades]);
title('Bin-blur vs Bin-blur-elements-removed');
subplot(223);
imshow([I_rotat, I_rotated_bin]);
title('I rotated vs I rot bin');
subplot(224);
imshow(I);
title('Original rotated + cropped');
pres_fig = gcf;
pres_fig.Units = 'normalized';
pres_fig.Position = [0 0 1 1];

function I_rem = remove_spades_numbers(I)
    [m,n] = size(I);
    counter_black = 0;
    counter_white = 0;
    change_from_white_to_black = 0;
    change_from_black_to_white = 0;
    total_pixels = m*n;
    
    white_indeces = find(I > 100 );
    for i=white_indeces(1):white_indeces(end)
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
    for i=2:total_pixels
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
            if counter_white > 0 && counter_white < 20
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

function [I_rotated, I_rotated_bin] = rotate_img(I_bin_rem,I_original)
    % Extract the Corner Indices
%% Extract the Corners
    [row,col] = find(I_bin_rem > 200);
    indices = [row,col];

    mins = min(indices);
    maxs = max(indices);

    minY = mins(1);
    minX = mins(2);
    maxY = maxs(1);
    maxX = maxs(2);

    [cornerX, cornerY] = find(I_bin_rem(minY,:) == 255);
    corner1 = [cornerY(1),minY];

    [cornerX, cornerY] = find(I_bin_rem(maxY,:) == 255);
    corner2 = [cornerY(1),maxY];

    [cornerX, cornerY] = find(I_bin_rem(:,minX) == 255);
    corner3 = [minX,cornerX(1)];

    [cornerX, cornerY] = find(I_bin_rem(:,maxX) == 255);
    corner4 = [maxX,cornerX(1)];

    thresholdDistanceX = 100;
    
    corners = [corner1; corner2; corner3; corner4];
    [~, x_min] = min(corners(:,1)); %lefmost corner x
    
    
    dist = zeros(4,1);
    for i=1:4
        dist(i) = vecnorm(corners(x_min,:) - corners(i,:),2,2);
    end
    %gets rid of largest distances
    [~,index_max]=max(dist);
    dist(index_max) = 0;
    [~,index_max]=max(dist);
    dist(index_max) = 0;
    %finds the smallest distance greater than zero(there's always a
    %distance which equals zero, corner - (same corner)
    [max_dist,index_max] = max(dist);
    
   %shortest side
    x = corners(index_max,1) - corners(x_min,1);
    y = corners(index_max,2) - corners(x_min,2);
    
    figure
    imshow(I_original);
    hold on
    text(corners(index_max,1),corners(index_max,2), num2str(index_max), 'FontSize', 30, 'Color', 'red');
    text(corners(x_min,1),corners(x_min,2), num2str(x_min), 'FontSize', 30, 'Color', 'blue');
    
    angle = rad2deg(atan2(y,x));
    %% Rotate Image
     % too close for rotation; distorts more slightly distorted images
    if (angle > 176 && angle < 184) || (angle > -4 && angle < 4)
        I_rotated = I_original;
        I_rotated_bin = I_bin_rem;
    else
        I_rotated = imrotate(I_original,angle,'crop');
        I_rotated_bin = imrotate(I_bin_rem,angle,'crop');
    end
end

function I_cropped = crop_img(I_bin_rotated, I_original_rotated)
    blackRowsXtop = 0;
    blackRowsXbottom = 0;
    blackRowsYleft = 0;
    blackRowsYright = 0;

    for i=1:size(I_bin_rotated,1)
        if sum(I_bin_rotated(i,:) == 255) > 0
            break
        else
            blackRowsXtop = blackRowsXtop+1;
        end
    end
    for i=blackRowsXtop:size(I_bin_rotated,1)
        if sum(I_bin_rotated(i,:) == 255) == 0
            blackRowsXbottom = blackRowsXbottom+1;
        end
    end

    for i=1:size(I_bin_rotated,2)
        if sum(I_bin_rotated(:,i) == 255) > 0
            break
        else
            blackRowsYleft = blackRowsYleft+1;
        end
    end
    for i=blackRowsYleft:size(I_bin_rotated,2)
        if sum(I_bin_rotated(:,i) == 255) == 0
            blackRowsYright = blackRowsYright+1;
        end
    end

    I_cropped = I_original_rotated(blackRowsXtop:size(I_bin_rotated,1)-blackRowsXbottom,...
        blackRowsYleft:size(I_bin_rotated,2)-blackRowsYright);
end

function I = check_if_upside_down(I)
    %% Check if the picture is upside down
    black_total_top_half =length(find(I(1:round(end/2),:) < 100));
    black_total_bottom_half =length(find(I(ceil(end/2):end,:) < 100));

    if black_total_bottom_half > black_total_top_half
        I = imrotate(I,180,'crop');
    end
end