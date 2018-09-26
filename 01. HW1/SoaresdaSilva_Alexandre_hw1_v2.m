%contrast enhancement
close all
clc, clear;
filter_int=@(n)1/(n^2)*ones(n);
L = 256;
I_original = imread('Testimage4.tif');
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
I_bin = binarize_img(I_blur,125);
%% remove spades / nums and rotate

I_bin_no_spades = remove_spades_numbers(I_bin);

figure
imshow([I_original, I_bin, I_bin_no_spades]);
[I_rotat, I_rotated_bin] = rotate_img(I_bin_no_spades,I_original);
I_cropped = crop_img(I_rotated_bin, I_rotat);
I = check_if_upside_down(I_cropped);
%% display figs

figure
%title('1. original VS 2. blur 5x5 VS 3. binarized VS 4. spades/no removed from binarized VS 5. sec derivative with contrast enhancement');
hold on
%imshow([I_original, I_blur, I_bin,I_rem, I_sec_deriv,I_rotat]);
subplot(221);
imshow([I_original, I_blur]);
subplot(222);
imshow([I_bin, I_bin_no_spades]);
subplot(223);
imshow([I_rotat, I_rotated_bin]);
subplot(223);
imshow([I_rotat, I_rotated_bin]);
subplot(224);
imshow([I]);
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
    
    [~,index_max]=max(dist);
    dist(index_max) = 0;
    [~,index_max]=max(dist);
    dist(index_max) = 0;
    [max_dist,index_max]=max(dist);
    
    
    %[~,index_min]=min(dist);
    %smallest side
    x = corners(index_max,1) - corners(x_min,1);
    y = corners(index_max,2) - corners(x_min,2);
    
    figure
    imshow(I_original);
    hold on
    text(corners(index_max,1),corners(index_max,2), num2str(index_max), 'FontSize', 30, 'Color', 'red');
    text(corners(x_min,1),corners(x_min,2), num2str(x_min), 'FontSize', 30, 'Color', 'blue');
    
    angle = rad2deg(atan2(y,x));
%     if angle < -120 || angle > 120
%         angle = angle  + 90;
%     end
    %imshow(imrotate(I_original,angle,'crop'))
        %% Rotate Image

    if (angle > 173 && angle < 183) || (angle > -5 && angle < 5)
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