%% Housekeeping
clc, clear, close all

run('VLFEATROOT/toolbox/vl_setup');
vl_version verbose

addpath('numbers');
load('digits_map.mat');
%% Capture Input Image
AVG_FILTER_SIZE = 3;
BIN = 0;


Ia_orig = cam_Capt;
Ia = extract_digits(Ia_orig,AVG_FILTER_SIZE,BIN);
Ia = Ia{1};
%% Predict and Match
thresh = 1;
[predictedNumberOptions, minScores, scores, normalized_scores,...
    no_matches, counts] =...
    matchAndPredict(digits_map, Ia, 1, thresh, AVG_FILTER_SIZE, BIN);
%% Functions
function I = cam_Capt

    webcamlist_conn = webcamlist;
    webcam_connected = webcam(webcamlist_conn{1});
%     webcam_connected.Resolution ='640x480';
    % webcam_connected.Resolution = '320x240';
    prv = preview(webcam_connected);
    prv.Parent.Visible = 'off';
%     prv.Parent.Parent.Parent.Position(1:2) = [0.6 .5];
    prv.Parent.Visible = 'on';
    prompt = 'Press enter to capture a frame';
    x = input(prompt);
    I = snapshot(webcam_connected);
    delete(prv);
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
    mean_g = mean2(I1);
    mean_g = mean_g/500;
    I_bin = ~imbinarize(uint8(I1),mean_g);
    for i=1:5
        I_bin = medfilt2(I_bin);
    end
    %% matlab nonuniform illumination
    I_cell = {I, I1, I_bin};
    I_props = regionprops(I_bin);

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
                
                
                if binarized
                    digit = imresize(digit_bin, [120 120]);
%                     se = strel('disk',2);
%                     digit = imdilate(digit,se);
%                     digit = (edge(digit,'Canny'));
                    digit = process_bin_num(digit,avg_filter_size);
                    digits{j} = digit;
                else
                    digit = imresize(imcrop(I, BoundingBox), [120 120]);
                    digits{j} = digit;
                end
                j = j + 1;
        end
    end    
end


function [digit, pre, pos] = I_crop_withBound(I_bin,boundingBox)
    digit = imcrop(I_bin, boundingBox);
    %             digits{i} = imcrop(I_binclosed,
    %             I_props(i).BoundingBox);close
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

function [predictedNumberOptions, minScores, scores,...
            normalized_scores, no_matches, counts] = ...
        matchAndPredict(templateImageMap, inputImage,... %function definition
            showMatches, thresh, AVG_FILTER_SIZE, BIN)
        
    digit = templateImageMap;
    Ia = inputImage;
    [fa,da] = vl_sift(im2single((Ia))) ;
    
    for j=1:length(digit)
        digit_no = digit(j-1);
        for i=1:length(digit_no)
            Ib_orig = digit_no{i};
            Ib = extract_digits(Ib_orig,AVG_FILTER_SIZE,BIN);
            Ib = Ib{1};
            [fb,db] = vl_sift(im2single((Ib))) ;
            [matches, score] = vl_ubcmatch(da,db,thresh);
            
            [matches, score] = nonUniqueReduction(matches,score);
            
            if(showMatches)
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
            end

            if(length(score)>1)
                median_score = median(score);
                indi_less_median = find(score<median_score);
                matches = matches(:,indi_less_median);
                score = score(indi_less_median);
            end

            scores(j,i) = sum(score);
            no_matches(j,i) = length(matches);
        end
    end

    for k=1:size(scores,2)
        normalized_scores(:,k) = abs((scores(:,k)-mean(scores(:,k)))./no_matches(:,k));
    end
    
    [minScores,predictedNumberOptions] = min(normalized_scores); 
    nans = isnan(minScores);
    minScores(nans) = inf;
    disp(['Predicted Numbers: ', num2str(predictedNumberOptions-1)])
    disp(['with scores: ', num2str(minScores)])
    
    counts = hist(predictedNumberOptions,1:10);
    measure = [];
    [count0,indiSort] = sort(counts,'descend');
    indiSort = indiSort-1;
    
    maxFreq1 = indiSort(1);
    maxFreq2 = indiSort(2);
    maxFreq3 = indiSort(3);
    
    countTemp = find(count0 ~= 0);
    count0 = count0(countTemp);
    
    if(all(diff(sort(count0(count0 ~= 0)))))
%         for k=0:9
%     %         frequency = counts(k);
%             indi = find(predictedNumberOptions-1 == k);
%             if(isempty(indi))
%                 measure(k+1) = inf;
%             else
%                 measure(k+1) = sum(minScores(indi))/length(indi);
%             end
%         end
%         [~,pI] = min(measure);

    disp(['Predicted Number: ', num2str(maxFreq1)])
    else
        [~,pI] = min(minScores);
        indi = find(predictedNumberOptions-1 == maxFreq1);
        score1 = sum(minScores(indi))/length(indi);

        indi = find(predictedNumberOptions-1 == maxFreq2);
        score2 = sum(minScores(indi))/length(indi);

        indi = find(predictedNumberOptions-1 == maxFreq3);
        score3 = sum(minScores(indi))/length(indi);


        disp(['Score: ', num2str(maxFreq1), ' - ' num2str(score1)])
        disp(['Score: ', num2str(maxFreq2), ' - ' num2str(score2)])
        disp(['Score: ', num2str(maxFreq3), ' - ' num2str(score3)])
    end
end

function [A3, B3] = nonUniqueReduction(matchings,scorings)
    Aa = matchings(1,:);
    Ab = matchings(2,:);
    Bc = scorings;

    n = length(Bc);
    while n ~= 0
        nn = 0;
        for ii = 1:n-1
            if Bc(ii) > Bc(ii+1)
                [Bc(ii+1),Bc(ii)] = deal(Bc(ii), Bc(ii+1));
                [Aa(ii+1),Aa(ii)] = deal(Aa(ii), Aa(ii+1));
                [Ab(ii+1),Ab(ii)] = deal(Ab(ii), Ab(ii+1));
                nn = ii;
            end
        end
        n = nn;
    end

    A2 = [Aa;Ab];
    B2 = Bc;

    Aa = A2(1,:);
    Ab = A2(2,:);
    Bc = B2;

    A3 = double.empty(2,0);
    B3 = double.empty(1,0);

    k = 1;
    for j = 1:length(Bc)
        aa = find(Aa == Aa(j));
        ab = find(Ab == Ab(j));
    
        idxs = [aa,ab];
        best = find(Bc == min(Bc(idxs)));
    
        if (isempty(find(A3 == Aa(best))) && isempty(find(A3 == Ab(best))))
            A3(1,k) = Aa(best);
            A3(2,k) = Ab(best);
            B3(k) = Bc(best);
            k = k + 1;
        end
    
    end
end

function I = process_bin_num(I_bin,avg_filter_size)
    LPF = avg_filt(avg_filter_size);
    I_bin = uint8(I_bin*255);
   
    % I1 = I;
%     I = conv2(I_bin,LPF,'same');
    I = I_bin;
    I = edge(I,'Canny');
    x = fspecial('motion',15,0);
%     y = fspecial('motion',25,90);
%     LPF = avg_filt(5);
%     
     I = conv2(I,x,'same');
end