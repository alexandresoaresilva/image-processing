clc, clear, close all

run('VLFEATROOT/toolbox/vl_setup');
vl_version verbose
addpath('new_ver');
load('Templates.mat');
load('descriptor_Templates.mat');
load('frame_Templates.mat');
load('image_Templates.mat');
load('digit.mat');

% imshow(I)
AVG_FILTER_SIZE = 3;

BIN = 0;
% compare_to = 6;
%% Pre-Process Image
Ia_orig = digit(1);
Ia_orig = Ia_orig{1};
% Ib_orig = image_Templates{compare_to};
% Ib_orig = cam_Capt;
% Ia_orig = imread('7_4.jpg');
% Ib_orig = imread('7_5.jpg');
Ia = extract_digits(Ia_orig,AVG_FILTER_SIZE,BIN);
% % Ib = Templates{8};
% Ib = extract_digits(Ib_orig,AVG_FILTER_SIZE,BIN);

Ia = Ia{1};
% Ib = Ib{1};

% imshow(Ia)

% 
% Ia = (edge(Ia,'Canny'));
% Ib = (edge(Ib,'Canny'));
% 
% se = strel('disk',2);
% Ia = imdilate(Ia,se);
% Ib = imdilate(Ib,se);

% Ia = imbinarize(imgaussfilt(double(Ia),2));
% Ib = imbinarize(imgaussfilt(double(Ib),2));

% figure(1)
% subplot(121)
% imagesc(Ia)
% colormap gray
% subplot(122)
% imagesc(Ib)
% colormap gray

%% SIFT 


[fa,da] = vl_sift(im2single((Ia))) ;
% [fb,db] = vl_sift(im2single((Ib))) ;

% figure(2)
% subplot(121)
% imshow(Ia)
% perm = randperm(size(fa,2)) ;
% sel = perm ;
% h1 = vl_plotframe(fa(:,sel)) ;
% h2 = vl_plotframe(fa(:,sel)) ;
% set(h1,'color','k','linewidth',3) ;
% set(h2,'color','y','linewidth',2) ;
% h3 = vl_plotsiftdescriptor(da(:,sel),fa(:,sel)) ;
% set(h3,'color','g') ;
% 
% figure(2)
% subplot(122)
% imshow(Ib);
% perm = randperm(size(fb,2)) ;
% sel = perm;
% h1 = vl_plotframe(fb(:,sel)) ;
% h2 = vl_plotframe(fb(:,sel)) ;
% set(h1,'color','k','linewidth',3) ;
% set(h2,'color','y','linewidth',2) ;
% h3 = vl_plotsiftdescriptor(db(:,sel),fb(:,sel)) ;
% set(h3,'color','g') ;

scores = [];
no_matches = [];
for i=1:length(Templates)
    Ib_orig = image_Templates{i};
    Ib = extract_digits(Ib_orig,AVG_FILTER_SIZE,BIN);
    Ib = Ib{1};
    fb = frame_Templates{i};
    db = descriptor_Templates{i};
%     db = Templates{i};
    [matches, score] = vl_ubcmatch(da,db,3);
    
    
    figure; clf ;
    imagesc(cat(2, Ia, Ib)) ;
    colormap gray
    axis image off ;

    xa = fa(1,matches(1,:)) ;
    xb = fb(1,matches(2,:)) + size(Ia,2) ;
    ya = fa(2,matches(1,:)) ;
    yb = fb(2,matches(2,:)) ;

    hold on ;
    h = line([xa ; xb], [ya ; yb]) ;
    set(h,'linewidth', 1, 'color', 'b') ;

    vl_plotframe(fa(:,matches(1,:)));
    fb(1,:) = fb(1,:) + size(Ia,2);
    vl_plotframe(fb(:,matches(2,:)));
    colormap gray
    axis image off ;

    
    [drop, perm] = sort(score, 'ascend') ;
    median_score = median(drop);
    drop(drop>median_score) = 0;
%     median_index = find()
%     perm = perm;
%     score  = score(perm) ;
    no_matches(i) = length(matches);
    scores(i) = sum(drop);
    
%     matches = matches(:, perm) ;
end

normalized_scores = abs((scores-mean(scores))./no_matches);
[~,pI] = min(normalized_scores); 
disp(['Predicted Number: ', num2str(pI-1)])


% 
% figure(3) ; clf ;
% imagesc(cat(2, Ia, Ib)) ;
% colormap gray
% axis image off ;
% 
% xa = fa(1,matches(1,:)) ;
% xb = fb(1,matches(2,:)) + size(Ia,2) ;
% ya = fa(2,matches(1,:)) ;
% yb = fb(2,matches(2,:)) ;
% 
% hold on ;
% h = line([xa ; xb], [ya ; yb]) ;
% set(h,'linewidth', 1, 'color', 'b') ;
% 
% vl_plotframe(fa(:,matches(1,:)));
% fb(1,:) = fb(1,:) + size(Ia,2);
% vl_plotframe(fb(:,matches(2,:)));
% colormap gray
% axis image off ;

% sum(scores)

%%
function I = cam_Capt

    % % webcam setup
    webcamlist_conn = webcamlist;
    webcam_connected = webcam(webcamlist_conn{1});
%     webcam_connected.Resolution ='640x480';
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

function digits = extract_digits(I,avg_filter_size,binarized)
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
    for i=1:5
        I_bin = medfilt2(I_bin);
    end
    %% matlab nonuniform illumination

    % se = strel('disk',STRUCTURING_ELEMENT_SIZE);
    % I_binclosed = imclose(I_bin,se);

    % I_cell = {I, I1, I_bin, I_binclosed};
    I_cell = {I, I1, I_bin};
    I_props = regionprops(I_bin);

    figure
%     imshow(I_bin);
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
                
                
                if binarized
                    digit = imresize(digit_bin, [120 120]);
                    digits{j} = digit;
                else
                    digit = imresize(imcrop(I, BoundingBox), [120 120]);
                    digits{j} = digit;
                end
                j = j + 1;
        end
    end    
end

function filt = avg_filt(n)
    filt = 1/(n^2)*ones(n);
end