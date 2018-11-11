%Database Creation
clc, clear, close all
image_Templates = cell(10,1);
descriptor_Templates = cell(10,1);
frame_Templates = cell(10,1);

i = 1;

while i < 11
    
    if i > 1
        figure(1)
        imshow(Ia,[]);
        perm = randperm(size(fa,2)) ;
        P = length(perm);
        sel = perm(1:P) ;
        h1 = vl_plotframe(fa(:,sel)) ;
        h2 = vl_plotframe(fa(:,sel)) ;
        set(h1,'color','k','linewidth',3) ;
        set(h2,'color','y','linewidth',2) ;
        h3 = vl_plotsiftdescriptor(da(:,sel),fa(:,sel)) ;
        set(h3,'color','g') ;
        figure(2)
        imshow(Ia,[]);
        x = input('Are you happy with the saved template? If not enter 0: ');
        if x == 0
            i = i - 1;
        else
            disp('Template saved');
        end
    end
    
    clc, clearvars -except image_Templates descriptor_Templates frame_Templates i, close all
    
    run('VLFEATROOT/toolbox/vl_setup');
    vl_version verbose
    
    AVG_FILTER_SIZE = 2;
    digit_size = 240;
    
    Ia_orig = cam_Capt;
    BIN = 0;
    Ia = extract_digits(Ia_orig,AVG_FILTER_SIZE,digit_size);
    
    Ia = Ia{1};
    
    [fa,da] = vl_sift(im2single((Ia))) ;
    
    image_Templates{i,1} = Ia;
    descriptor_Templates{i,1} = da;
    frame_Templates{i,1} = fa;
    
    i = i + 1;
    
end

save('image_Templates.mat', 'image_Templates');
save('descriptor_Templates.mat', 'descriptor_Templates');
save('frame_Templates.mat', 'frame_Templates');


%%
function I = cam_Capt

% % webcam setup
webcamlist_conn = webcamlist;
webcam_connected = webcam(webcamlist_conn{1});
% webcam_connected.Resolution ='640x480';
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

function [digits, digits_bin] = extract_digits(I,avg_filter_size,digit_side)
    STRUCTURING_ELEMENT_SIZE = 2;
    
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
    I_bin = ~imbinarize(uint8(I1),.5);
    for i=1:20
        I_bin = medfilt2(I_bin);
    end
    %% matlab nonuniform illumination
    se = strel('disk',STRUCTURING_ELEMENT_SIZE);
    I_bin = imclose(I_bin,se);

    I_cell = {I, I1, I_bin};
    I_props = regionprops(I_bin);

%     figure
%     imshow(I_bin);
%     hold on;
    digits = cellmat(0);
    j = 1;
    for i=1:length(I_props)
        if I_props(i).Area > (M*N / 330)
%                 rect = rectangle('Position',I_props(i).BoundingBox,...
%                     'EdgeColor','r','LineWidth',3);
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

function [digit, pre, pos] = I_crop_withBound(I_bin,boundingBox)
    digit = imcrop(I_bin, boundingBox);
    
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

function filt = avg_filt(n)
    filt = 1/(n^2)*ones(n);
end