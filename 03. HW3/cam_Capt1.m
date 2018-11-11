clc, clear, close all
addpath('VLFEATROOT')
run('VLFEATROOT/toolbox/vl_setup');
vl_version verbose
%for the closing operation
% SELECT_I = 8;
% 
% AVG_FILTER_SIZE = 3;
% [file_names, file_names_char] = save_file_names_in_folder(pwd,'jpg');

% imshow(I)
AVG_FILTER_SIZE = 3;
DIGIT_SIZE = 120;
% BIN = 0
%% Pre-Process Image
% Ia_orig = cam_Capt; %capturing image
% Ib_orig = cam_Capt;
a = load('image_Templates.mat');
Ia_orig = a.image_Templates{7}
% Ia_orig = imread('6.jpg');
Ib_orig = imread('7.jpg');

% [Ia, Ia_bin] = extract_digits(Ia_orig,AVG_FILTER_SIZE,DIGIT_SIZE);
Ia_bin = ~imbinarize(Ia_orig,.35);
Ia = Ia_bin;
[Ib, Ib_bin] = extract_digits(Ib_orig,AVG_FILTER_SIZE,DIGIT_SIZE);

% Ia = Ia{1};
% Ib = Ib{1};
% Ia = Ia_bin{1};
Ib = Ia_bin;
% Ib = Ib_bin{1};

imshow(Ia)
% Ia = (edge(Ia,'Canny'));
% Ib = (edge(Ib,'Canny'));
% 
% se = strel('disk',1);
% Ia = imdilate(Ia,se);
% Ib = imdilate(Ib,se);

% Ia = imbinarize(imgaussfilt(double(Ia),2));
% Ib = imbinarize(imgaussfilt(double(Ib),2));

figure(1)
subplot(121)
imagesc(Ia)
colormap gray
subplot(122)
imagesc(Ib)
colormap gray
%% SIFT 
[fa,da] = vl_sift(im2single(Ia)) ;
[fb,db] = vl_sift(im2single(Ib)) ;
a = load('Templates.mat')
da = a.Templates{2};

figure(2)
subplot(121)
imshow(Ia)
perm_indeces = randperm(size(fa,2)) ;

sel = perm_indeces ;
h1 = vl_plotframe(fa(:,sel)) ;
h2 = vl_plotframe(fa(:,sel)) ;
set(h1,'color','k','linewidth',3) ;
set(h2,'color','y','linewidth',2) ;

h3 = vl_plotsiftdescriptor(da(:,sel),fa(:,sel)) ;
set(h3,'color','g') ;

figure(2)
subplot(122)
imshow(Ib);
perm_indeces = randperm(size(fb,2)) ;
sel = perm_indeces;
h1 = vl_plotframe(fb(:,sel)) ;
h2 = vl_plotframe(fb(:,sel)) ;
set(h1,'color','k','linewidth',3) ;
set(h2,'color','y','linewidth',2) ;
h3 = vl_plotsiftdescriptor(db(:,sel),fb(:,sel)) ;
set(h3,'color','g') ;

[matches, scores] = vl_ubcmatch(da,db);
[drop, perm_indeces] = sort(scores, 'ascend') ;
% perm_indeces = perm_indeces;
matches = matches(:, perm_indeces) ; %sorts matches
scores  = scores(perm_indeces) ;  %sorts scores

figure(3) ; clf ;
imagesc(cat(2, Ia, Ib)) ;
colormap gray
axis image off ;

% angles_diff = rad2deg(abs(fa(4,matches(1,:))-fb(4,matches(2,:))));
% median_angles_diff = median(angles_diff)
% j=1; %index for x,y
% selected_matches = zeros(3,1);
% for i=1:length(angles_diff)
%     if angles_diff(i) < (median_angles_diff+45)
%         selected_matches(:,j) = matches(:,i);
%         xa(j) = fa(1,matches(1,i)) % 1st row of matches are indeces of fa ;
%         xb(j) = fb(1,matches(2,i)) + size(Ia,2) ; % 2nd row of matches are indeces of fb ;
%         ya(j) = fa(2,matches(1,i)) ;
%         yb(j) = fb(2,matches(2,i)) ;
%         j = j + 1;
%     end 
% end
xa = fa(1,matches(1,:)) % 1st row of matches are indeces of fa ;
xb = fb(1,matches(2,:)) + size(Ia,2) ; % 2nd row of matches are indeces of fb ;
ya = fa(2,matches(1,:)) ;
yb = fb(2,matches(2,:)) ;

pairsA = [xa' ya'];
pairsB = [xb' yb']; %if it is unique and has the same magnitude & angle

BOUNDARY = 5;
seen_before = 0;
selected_matches = zeros(2,1);
newscores=0;
close_indeces = 0;
k = 1;
scores_find = 0;
i = 1;
index_found_matches = 0;
while(i)
    lower_b = pairsB(i,:)-BOUNDARY;
    upper_b = pairsB(i,:)+BOUNDARY ;
    first_index = 0;
    seen_before = 0;
    close_indeces = 0;
    k = 1;
    %run through the whole list to find similar x, y's
    for j=1:length(pairsB) 
        if ( pairsB(j,1) > lower_b(1) && pairsB(j,1) < upper_b(1))
            if (pairsB(j,2) > lower_b(2) && pairsB(j,2) < upper_b(2))
                if ~seen_before %NOT seen before
                    seen_before = 1;
                    first_index = j;
                    index_found_matches = index_found_matches + 1;
                end
                close_indeces(k) = j;
                k = k + 1;
            end
        end %if pairsB
    end
    close_scores = scores(close_indeces);
    [smallest_score, index_min_close] = min(close_scores(close_scores>0));
    index_min = find(scores == smallest_score);
    %erases other close pairs and saves the minimum score
%     pair_min = pairsB(index_min,:);
%     pairsB(close_indeces,:) = 0;
%     pairsB(index_min,:) = pair_min;
    
    scores(close_indeces) = 0;
    scores(index_min) = smallest_score;
    
    selected_matches(1:2,index_found_matches) = matches(:,index_min); %saves the i_th selection
%     selected_matches(3,k) = j; %saves the corresponding index of the match
    clear close_indeces;
%      pairsB(:,close_indeces) = zeros(2,length(close_indeces));
    seen_before = 0;
    
    i = i + 1;
    if k == length(pairsB) || i == length(pairsB)
        break;
    end
end
hold on ;

h = line([xa ; xb], [ya ; yb]) ;
set(h,'linewidth', 1, 'color', 'b') ;

vl_plotframe(fa(:,selected_matches(1,:))) ;
fb(1,:) = fb(1,:) + size(Ia,2) ;
vl_plotframe(fb(:,selected_matches(2,:))) ;
colormap gray
axis image off ;

sum(scores);
%%
function remove_bkgrnd_matches(I,I_bin,f, d)

        remove_bkgrnd_matches

end
function I = cam_Capt

    % % webcam setup
    webcamlist_conn = webcamlist;
    webcam_connected = webcam(webcamlist_conn{1});
    webcam_connected.Resolution ='640x480';
    % webcam_connected.Resolution = '320x240';
    prv = preview(webcam_connected);
    prv.Parent.Visible = 'off';
    prv.Parent.Parent.Parent.Position(1:2) = [0.6 .5];
    prv.Parent.Visible = 'on';
    prompt = 'Press enter to capture a frame';
    x = input(prompt);
    I = snapshot(webcam_connected);
    delete(prv);
end
function [digit, pre, pos] = I_crop_withBound(I_bin,boundingBox)
    digit = imcrop(I_bin, boundingBox);
    %             digits{i} = imcrop(I_binclosed, I_props(i).BoundingBox);
    [m,n] = size(digit);
    [max_x,i] = max([m,n]);            
    [digit, pre, pos] = pad_digit_img(digit ,ceil(max_x*1.05));

    
    function [I, pre, pos] = pad_digit_img(I,square_side)
        [x,y] = size(I);
        [x_pre, x_pos] = padding_pre_post_calc(x,square_side);
        [y_pre, y_pos] = padding_pre_post_calc(y,square_side);

        pre = [x_pre y_pre];
        pos = [x_pos y_pos];
        I = padarray(I, pre,'pre');
        I = padarray(I, pos,'pos');
        
        function [n_pre, n_pos] = padding_pre_post_calc(n,dim_size)
            n_is_not_even = rem(dim_size-n,2);
            n_pre = (dim_size-n)/2;
    %         n_pre_is_not_even = rem(n_pre,2);

            if n_is_not_even && n_pre %n_pre not zero
                n_pre = n_pre - 0.5;
                n_pos = n_pre + 1;
            else
                n_pos = n_pre;
            end
            n_pos = double(n_pos);
            n_pre = double(n_pre);
        end
    end
end

function [digits, digits_bin] = extract_digits(I,avg_filter_size,digit_side)
    if length(size(I)) > 2
        I = rgb2gray(I);
    end
    %% procesing
    % I = imresize(I,[240 320]);
    [M, N] = size(I);
    LPF = avg_filt(avg_filter_size);
    % I1 = I;
    I1 = conv2(I,LPF,'valid');
    I_bin = I1;
    I_bin = ~imbinarize(uint8(I1),.35);
    for i=1:20
        I_bin = medfilt2(I_bin);
    end
    %% matlab nonuniform illumination

    % se = strel('disk',STRUCTURING_ELEMENT_SIZE);
    % I_binclosed = imclose(I_bin,se);

    % I_cell = {I, I1, I_bin, I_binclosed};
    I_cell = {I, I1, I_bin};
    I_props = regionprops(I_bin);

    figure
    imshow(I_bin);
    hold on;
    digits = cellmat(0);
    j = 1;
    for i=1:length(I_props)
        if I_props(i).Area > (M*N / 330)
                rect = rectangle('Position',I_props(i).BoundingBox,...
                    'EdgeColor','r','LineWidth',3);
                [digit_bin, pre, pos] = I_crop_withBound(I_bin,I_props(i).BoundingBox);
                
                x_origin = I_props(i).BoundingBox(1)-pre(2);
                y_origin = I_props(i).BoundingBox(2)-pre(1);
                
                x_width = I_props(i).BoundingBox(3) + pre(2)+pos(2);
                y_width = I_props(i).BoundingBox(4) + pre(1)+pos(1);
                BoundingBox = [x_origin,y_origin,...
                    x_width, y_width];
                
                digitBin = imresize(digit_bin, [digit_side digit_side]);
                digits_bin{j} = digitBin;

                digit = imresize(imcrop(I, BoundingBox), [digit_side digit_side]);
                digits{j} = digit;
                
                j = j + 1;
        end
    end    
end

function filt = avg_filt(n)
    filt = 1/(n^2)*ones(n);
end