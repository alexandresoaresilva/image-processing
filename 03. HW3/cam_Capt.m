clc, clear, close all
% addpath('VLFEATROOT')
run('VLFEATROOT/toolbox/vl_setup');
vl_version verbose
%for the closing operation
SELECT_I = 8;

AVG_FILTER_SIZE = 3;
[file_names, file_names_char] = save_file_names_in_folder(pwd,'jpg');

% webcam setup
% webcamlist_conn = webcamlist;
% webcam_connected = webcam(webcamlist_conn{1});
% webcam_connected.Resolution ='640x480';
% 
% webcam_connected.Resolution = '320x240';
% prv = preview(webcam_connected);
% prv.Parent.Visible = 'off';
% prv.Parent.Parent.Parent.Position(1:2) = [0.6 .5];
% prv.Parent.Visible = 'on';
% prompt = 'Press enter to capture a frame';
% x = input(prompt);
% I = snapshot(webcam_connected);
% delete(prv);
%% load file

I = imread('fuzzy3_2.jpg');
digits = extract_digits(I,AVG_FILTER_SIZE);

I2 = imread('IMG_20181102_100540197.jpg');
digits2 = extract_digits(I2,AVG_FILTER_SIZE);

figure;
sub_plots(digits);
%% sift code
% run('VLFEATROOT/toolbox/vl_setup');
% vl_version verbose

[f,d] = vl_sift(single(digits{1}));

perm = randperm(size(f,2)) ;
P = length(perm);
sel = perm(1:P);
h1 = vl_plotframe(f(:,sel)) ;
h2 = vl_plotframe(f(:,sel)) ;
set(h1,'color','k','linewidth',2) ;
set(h2,'color','y','linewidth',1) ;
h3 = vl_plotsiftdescriptor(d(:,sel),f(:,sel)) ;
set(h3,'color','g') ;


% I2 = imgaussfilt(I2,20);
% I2 = wiener2(I2,[10 10]);
% I2 = single(edge(I2,'Canny'));
% 
% I2 = imclose(I2,se);

% subplot(1,2,2)
% imshow(I2,[]);

[f1,d1] = vl_sift(single(digits2{7}));
perm = randperm(size(f1,2)) ;
P = length(perm);
sel = perm(1:P);
h1 = vl_plotframe(f1(:,sel)) ;
h2 = vl_plotframe(f1(:,sel)) ;
set(h1,'color','k','linewidth',2) ;
set(h2,'color','y','linewidth',1) ;
h3 = vl_plotsiftdescriptor(d1(:,sel),f1(:,sel)) ;
set(h3,'color','g') ;

[matches, scores] = vl_ubcmatch(d, d1);
%% functions

%digits is a cell array
function digits = extract_digits(I,avg_filter_size)
    if length(size(I)) > 2
        I = rgb2gray(I);
    end
    %% procesing
    % I = imresize(I,[240 320]);
    [M, N] = size(I);
    LPF = avg_filt(avg_filter_size);
    % I1 = I;
    I1 = conv2(I,LPF,'valid');

    I_bin = ~imbinarize(uint8(I1),.4);
    for i=1:20
        I_bin = medfilt2(I_bin);
    end
    %% matlab nonuniform illumination

    % se = strel('disk',STRUCTURING_ELEMENT_SIZE);
    % I_binclosed = imclose(I_bin,se);

    % I_cell = {I, I1, I_bin, I_binclosed};
    I_cell = {I, I1, I_bin};
    I_props = regionprops(I_bin);

%     figure
%     imshow(I_bin);
    hold on;
    digits = cellmat(0);
    j = 1;
    for i=1:length(I_props)
        if I_props(i).Area > (M*N / 330)
            if I_props(i).BoundingBox(4) > N/8 %n is height
                rect = rectangle('Position',I_props(i).BoundingBox,'EdgeColor','r','LineWidth',3);
                digit = I_crop_withBound(I_bin,I_props(i).BoundingBox);
                digits{j} = I_close(digit,150);
                j = j + 1;
            end
        end
    end    
end

function digit = I_close(digit, resize_side)
    h = histogram(digit);
    h = h.Values(2); %white values
    a = length(digit);
    avg_Width = floor(h/a);
    structuring_element_size = floor(avg_Width/4);
    se = strel('disk',structuring_element_size);
    digit = imclose(digit,se);
    digit = imresize(digit,[resize_side resize_side]);
end

function digit = I_crop_withBound(I_bin,boundingBox)
    digit = imcrop(I_bin, boundingBox);
    %             digits{i} = imcrop(I_binclosed, I_props(i).BoundingBox);
    [m,n] = size(digit);
    [max_x,i] = max([m,n]);            
    digit = pad_digit_img(digit ,ceil(max_x*1.2));

    
    function I = pad_digit_img(I,square_side)
        [x,y] = size(I);
        [x_pre, x_pos] = padding_pre_post_calc(x,square_side);
        [y_pre, y_pos] = padding_pre_post_calc(y,square_side);

        I = padarray(I,[x_pre y_pre],'pre');
        I = padarray(I,[x_pos y_pos],'pos');
        
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

function sub_plots(I_cell)
    n = length(I_cell);
    no_of_plots_per_row = ceil(n/2);
%     figure;
    for i=1:n
        subplot(2,no_of_plots_per_row,i)
        imshow(I_cell{i},[]);
    end
end

function [I2, I3, I4] = matlab_nonuni_illu_remov(I,disk_size, bin_thresh)
    se = strel('disk',disk_size);
    background = imopen(I,se);
%     LPF = avg_filt(50);
%     background = conv2(I,LPF,'same');
    %imshow(background);
    I2 = I - background;
    %imshow(I2,[]);
    I3 = imadjust(I2);
    %imshow(I3,[]);
    bw = imbinarize(I3);
    I4 = bwareaopen(bw,bin_thresh);
end

function filt = avg_filt(n)
    filt = 1/(n^2)*ones(n);
end