clc, clear, close all
addpath('vlfeat-0.9.21')

run('vlfeat-0.9.21/toolbox/vl_setup');
addpath('numbers'); %where images/mat files are stored are stored
load('digits_map.mat');
%% paramteters
AVG_FILTER_SIZE = 3; %for extracting digits
DIGIT_SIZE = 120; %image size
MOTION = 24;
MATCH_THRESHOLD = 2.2;
BINARIZE = 1;
SHOW_MATCHES = 0;
DIGIT_SELECTED = 9;
% most_freq_num

% AVG_FILTER_SIZE = 3; %for extracting digits
% DIGIT_SIZE = 120; %image size
% MOTION = 24;
% MATCH_THRESHOLD = 2.2;
% BINARIZE = 1;
% SHOW_MATCHES = 1;
% DIGIT_SELECTED = 9;
% %% paramteters
% numbers are predicted 6/10 with these parameters for alex's (
% using dig{8}, dig{5}, dig{2}})
% dig{1}; 1/10
% dig{3}; 4/10
% dig{4}; 5/10
% dig{6}; 1/10
% dig{7}; 4/10


Ia_orig = cam_capt; 
%alex version
[most_freq_num, ~, ~, total_matches] =...
    match_and_predict(Ia_orig, digits_map,SHOW_MATCHES, AVG_FILTER_SIZE,...
                      MATCH_THRESHOLD, DIGIT_SIZE, BINARIZE, MOTION)


% AVG_FILTER_SIZE = 3; %for extracting digits
% DIGIT_SIZE = 120; %image size
% MOTION = 24;
% MATCH_THRESHOLD = 2.2;
% BINARIZE = 1;
% SHOW_MATCHES = 0;
% DIGIT_SELECTED = 9;
%% Loading and pre-processing images
%image A & B
function run_test_digits()
    accumulate = [];
    avg_Mohak = [];
    avg_Alex = [];
% for i=1:10
%     DIGIT_SELECTED = i-1;
%      dig = digits_map(DIGIT_SELECTED);
%     Ia_orig = dig{8};   
    Ia_orig = cam_capt; 
    
    disp('                                           ');
    disp('===========================================');
    disp(['RUN & number: ', num2str(DIGIT_SELECTED)]);
    disp('-------------------------------------------');
%     for j=1:8
        
        [predictedNumberOptions, minScores, scores, normalized_scores,...
            no_matches, counts] = matchAndPredict(digits_map, Ia_orig, 0,...
            MATCH_THRESHOLD, AVG_FILTER_SIZE, DIGIT_SIZE, 1);
        
        %alex version
        [most_freq_num, scores_sum,no_matches, total_matches] =...
            match_and_predict(Ia_orig, digits_map,SHOW_MATCHES, AVG_FILTER_SIZE,...
                              MATCH_THRESHOLD, DIGIT_SIZE, BINARIZE, MOTION);
        disp('-------------------------------------------');
        disp(['number: ', num2str(DIGIT_SELECTED)]);
%         disp(['PREDICTED Mohak: ', num2str(predictedNumberOptions-1)]); 
        disp(['PREDICTED Alex: ', num2str(most_freq_num)]);
%         if most_freq_num == (i-1)
%             avg_Alex(i) = 1;
%         else
%              avg_Alex(i) = 0;
%         end
% %         
%         if ((predictedNumberOptions-1) == i)
%             accumulate(j,1) = 1;
%         else
%             accumulate(j,1) = 0;
%         end
%         if (most_freq_num == i)
%             accumulate(j,2) =  1;
%         else
%             accumulate(j,2) = 0;
%         end
% %     end
%     avg_Mohak(i+1) = mean(accumulate(:,1))*100;
%     avg_Alex(i+1) = mean(accumulate(:,2))*100;
%     disp('-------------------------------------------');
%     disp(['number: ', num2str(DIGIT_SELECTED)]);
%     disp(['acc. Mohak after 8 runs: ', num2str(avg_Mohak(i+1)),' %']); 
%     disp(['acc. Alex after 8 runs: ', num2str(avg_Alex(i+1)),' %']);
% end
disp(['accuracy: ', num2str(mean(avg_Alex)*100),' %']);
% avg_Mohak = mean(avg_Mohak)*100;
% avg_Alex = mean(avg_Alex)*100;
% disp('-------------------------------------------');
% disp(['number: ', num2str(DIGIT_SELECTED)]);
% disp(['acc. Mohak total: ', num2str(avg_Mohak),' %']); 
% disp(['acc. Alex total: ', num2str(avg_Alex),' %']);
end

function [predictedNumberOptions, minScores, scores, normalized_scores,...
    no_matches, counts] = ...
    matchAndPredict(templateImageMap, Ia, showMatches, SIFT_match_thresh,...
                    avg_filter_size,digit_size, binarized)
        
    digit = templateImageMap;
    [fa,da] = vl_sift(im2single((Ia))) ;
    
%     for j=1:length(digit)
    for j=1:10
        digit_no = digit(j-1);
        
        for i=1:length(digit_no)
            Ib_orig = digit_no{i};
%             [~, Ibin] = extract_digits(Ib_orig, avg_filter_size, digit_size);
%             Ib = process_bin_num(Ib{1},2);
            [Ib, Ibin] = extract_digits(Ib_orig, avg_filter_size, digit_size);
            
            if binarized
              	Ib = process_bin_num(Ibin{1},avg_filter_size, 15);
            else
                Ib = Ib{1};
            end
            
%             [fb,db] = vl_sift(im2single(Ib),'edgethresh', 3.5,'PeakThresh',10) ;
            [fb,db] = vl_sift(im2single(Ib)) ;
            [matches, score] = vl_ubcmatch(da,db,SIFT_match_thresh);
            [matches, score] = clean_features(matches, score, fa, fb,size(Ib,2));
            if(showMatches)
                figure; clf ;
                plot_SIFT_lines(Ia, Ib, fa, fb, matches);
            end
            %colecting the sum of scores to find the smallest among the
            %most frequent 
%             scores(j,i) = sum(score);
            no_matches(j,i) = length(matches);
            norm_scores = normalize_scores(score);
            if no_matches(j,i) > 0
                scores(j,i) = sum(norm_scores)/no_matches(j,i);
            else
                scores(j,i) = 1;
            end
        end
    end
    
%     for k=1:size(scores,2)
%         normalized_scores(:,k) = abs((scores(:,k)-mean(scores(:,k)))./no_matches(:,k));
%     end
    normalized_scores = scores;
    [minScores, predictedNumberOptions] = min(normalized_scores); 
%     nans = isnan(minScores);
%     minScores(nans) = inf;
    disp(['Predicted Numbers: ', num2str(predictedNumberOptions-1)]);
    disp(['with scores: ', num2str(minScores)]);
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

function [most_freq_num, scores_sum,no_matches, total_matches] =...
    match_and_predict(I, digits_map,show_matches, avg_filter_size,...
                      match_thresh, digit_size, binarize, motion_size)
    
    [Ia, Ia_bin] = extract_digits(I,avg_filter_size,digit_size);
    if binarize
       Ia = process_bin_num(Ia_bin{1},avg_filter_size, motion_size);
    else
       Ia = Ia{1};
    end
    [fa,da] = vl_sift(single(Ia)) ;
    scores_sum  = []; 
    no_matches = [];
    total_matches = [];
    %looks for the digit match
%     for i=1:length(digits_map)
    for i=1:10
        test_dig = digits_map(i-1);
        if show_matches
            figure;
        end
        for j=1:length(test_dig)
           [Ib, Ib_bin] =...
               extract_digits(test_dig{j}, avg_filter_size, digit_size);
           if binarize %process it in binary with motion blur
               Ib = process_bin_num(Ib_bin{1},avg_filter_size, motion_size);
           else
               Ib = Ib{1};
           end
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %% SIFT
            [fb,db] = vl_sift(single(Ib));
            [matches, scores] = vl_ubcmatch(da, db, match_thresh);
            [matches, scores] = clean_features(matches, scores, fa, fb,size(Ib,2));
            
            if show_matches
                subplot(2,4,j);
                plot_SIFT_lines(Ia, Ib, fa, fb, matches);                
            end
           
            norm_scores = normalize_scores(scores);
            no_matches(i,j) = length(matches);
            if no_matches(i,j) > 0
                scores_sum (i,j) = sum(norm_scores)/no_matches(i,j);
            else
                scores_sum (i,j) = 1;
            end
            
        end
        total_matches(i) = sum(no_matches(i,:));
    end
    [~,k] = min(sum(scores_sum,2));
    most_freq_num = k-1;
%      most_freq_num = find_match(scores_sum, total_matches);
%     most_freq_num = find_most_frequent_digit(scores_sum);
end

function [matches, scores] = clean_features(matches, scores, fa, fb,offset)
    [matches, scores] = remove_inf_and_nans(matches, scores);
    [matches, scores] = unique_matches_scores(fb,0, matches, scores);    
    [matches, scores] = unique_matches_scores(fa,1,matches, scores);
    [matches, scores] = remove_senseless_scores(matches, scores, fa, fb,offset);
    [matches, scores] = remove_outliers(matches, scores, 0);
%     [matches, scores] = keep_only_three_best_features(matches, scores);
end

function [matches,scores] = keep_only_three_best_features(matches, scores)
    if length(scores) > 3
        [score1, i1] = min(scores);
        scores(i1) = 0;
        [score2, i2] = min(scores(scores>0));
        scores(i2) = 0;
        [score3, i3] = min(scores(scores>0));
        scores = [score1 score2 score3];
        matches = matches(:,[i1 i2 i3]);
    else
        scores = scores ;
        matches = matches;      
    end
end

function i_scr = find_match(scores_sum, total_matches)
    sum_columns_scores = sum(scores_sum,2);
    [~, i_scr] = min(sum_columns_scores);
    [~, i_total]= max(total_matches);
    
    for k=1:length(total_matches)
        if i_total ~= i_scr
            total_matches(i_total) = 0;
            sum_columns_scores(i_scr) = 10;
        else
            break;
        end
        [~, i_total]= max(total_matches);    
        [~, i_scr] = min(sum_columns_scores);
    end
    i_scr = i_scr - 1;
end

function most_freq_num = find_most_frequent_digit(scores_sum)
    [min_scrs,i_min_scr] = min(scores_sum);
    digits_found_freq = i_min_scr-1;
    tab = tabulate(digits_found_freq);
    [~,index] = max(tab(:,2)); %gets number with largest count
    most_freq_num = tab(index,1);
end

function norm_scores = normalize_scores(scores)
     L = length(scores);
    if  L > 1
        min_scr = min(scores);
        max_scr = max(scores);
        norm_scores = (scores - min_scr)./(max_scr - min_scr);
    else
        norm_scores = scores/scores;
    end
end

function [matches, scores] = remove_inf_and_nans(matches, scores)
    i_nan = find(isnan(scores));
    if length(i_nan) > 0
        i_NOT_nan = find(~isnan(scores));
        scores = scores(i_NOT_nan);
        matches = matches(i_NOT_nan);
    end
    
    i_inf = find(isinf(scores));
    if length(i_inf) > 0
        i_NOT_inf = find(~i_inf(scores));
        scores = scores(i_NOT_inf);
        matches = matches(i_NOT_inf);
    end
end

function [matches, scores] = remove_outliers(matches, scores, no_std)
    
    if length(scores) > 4
        med_scr = median(scores);
        std_scr = std(scores);
        i = find(scores < (med_scr + std_scr*no_std ) );
        scores = scores(i);
        matches = matches(:,i);
    else
        scores = scores;
        matches = matches;
    end
end

function [matches,scores] = unique_matches_scores(fb_or_fa, is_fa,...
                                                matches,scores)
    index = 2;
    if is_fa
        index = 1;
    end
%     xa = fa(1,matches(1,:)) ;
    xb = fb_or_fa(1,matches(index ,:)) ;
%     ya = fa(2,matches(1,:)) ;
    yb = fb_or_fa(2,matches(index,:)) ;
    [xb_un, ixa, ixc] = unique(round(xb));
    [yb_un, iya, iyc] = unique(round(yb));
    
    if length(xb_un) >= length(yb_un)
       matches=  matches(:,iya);
       scores=  scores(:,iya);
    else
       matches=  matches(:,ixa);
       scores=  scores(:,ixa);
    end
end

function [matches, scores] = remove_senseless_scores(matches, scores, fa, fb, offset)
    L = length(scores);
    
    if L == 1
        scores = scores;
        matches = matches;
    elseif L == 2
        [~,ind] = min(scores);
        scores = scores(ind);
        matches = matches(:,ind);
    else
        xa = fa(1,matches(1,:)) ;
        xb = fb(1,matches(2,:))+offset;
        ya = fa(2,matches(1,:)) ;
        yb = fb(2,matches(2,:)) ;

        slopes = (ya-yb)./(xa-xb);

        indeces_pos = find(slopes >= 0);
        indeces_neg = find(slopes < 0);

        if length(indeces_pos) > length(indeces_neg)
            matches = matches(:, indeces_pos);
            scores = scores(indeces_pos);
        elseif length(indeces_pos) < length(indeces_neg)
            matches = matches(:, indeces_neg);
            scores = scores(indeces_neg);
        else
            matches = matches;
            scores = scores;
        end
    end

end

function [unique_maches, unique_scores] = nonUniques(matchings,scorings)
    Aa = matchings(1,:);
    Ab = matchings(2,:);
    Bc = scorings;
    if unique(scorings) == 0 %it's the same picture
        unique_maches = matchings;
        unique_scores = scorings;
    else
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

        unique_maches = double.empty(2,0);
        unique_scores = double.empty(1,0);
        k = 1;
        for j = 1:length(B2)           %Reduce to best unique matches
            P = A2(1,j);
            if(isempty(find(unique_maches == P)))
                Q = A2(2,j);
                idxOfUniqueMinScore = find(B2 == min(B2(round((find(A2 == Q).')/2))));
                idxOfUniqueMinMatch = A2(2,idxOfUniqueMinScore(1));
                if(isempty(find(unique_maches == Q)))
                    unique_maches(1,k) = A2(1,j);
                    unique_maches(2,k) = A2(2,idxOfUniqueMinScore);
                    unique_scores(k) = B2(idxOfUniqueMinScore(1));
                    k = k + 1;
                end

            end

        end
    end
end

function predict_no(i, scores_collected, no_matches)

    normalized_scores = abs(scores_collected-mean(scores_collected));
    [~,pI] = min(normalized_scores);
    normalized_scores(pI) = 0;
    [~,pI_2] = min(normalized_scores(normalized_scores>0));
    disp(['Real: ',num2str(i), ' Predicted Number: ', num2str(pI-1)])
    disp(['Next best match: ', num2str(pI_2-1)])

end

function [f,d] = run_SIFT(I, x_offset)
    % A frame is a disk of center f(1:2), scale f(3) and orientation f(4)
    [f,d] = vl_sift(single(I)) ;
    perm = randperm(size(f,2)) ;
    sel = perm;
    %two circles because one is larger thant the other and thus, becomes
    %   a little contour (plotting circles with angles
    % the offset is present for the second picture
    f_for_plotting = f;
    f_for_plotting(1,:) = f_for_plotting(1,:) + x_offset;
    
%     h1 = vl_plotframe(f_for_plotting(:,sel)) ;
%     h2 = vl_plotframe(f_for_plotting(:,sel)) ;
%     set(h1,'color','k','linewidth',3) ;
%     set(h2,'color','y','linewidth',2) ;
    
    %plots square frame 4x4 with gradients
%     h3 = vl_plotsiftdescriptor(d(:,sel),f_for_plotting(:,sel)) ;
%     set(h3,'color','g') ;
end

function plot_SIFT_lines(Ia, Ib, fa, fb, matches)
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

function I = process_bin_num(I_bin,avg_filter_size, motion_blur_amount)
    
    I = uint8(I_bin*255);
    if motion_blur_amount
        I = edge(I,'Canny');
        x = fspecial('motion',motion_blur_amount,0);
        I = conv2(I,x,'valid');        
    else
        LPF = avg_filt(avg_filter_size);
        I = conv2(I_bin,LPF,'valid');
    end
end

function I = process_I(I)
%     LPF = avg_filt(2);
%     I = uint8(I*255);
   
    % I1 = I;
%     I = conv2(I,LPF,'same');
%     I = edge(I,'Canny');
    x = fspecial('motion',15,0);
%     y = fspecial('motion',25,90);
%     LPF = avg_filt(5);
%     
     I = conv2(I,x,'same');
end

function filt = avg_filt(n)
    filt = 1/(n^2)*ones(n);
end

function [digits, digits_bin] = extract_digits(I, avg_filter_size, digit_side)
    STRUCTURING_ELEMENT_SIZE = 3;
    
    if length(size(I)) > 2
        I = rgb2gray(I);
    end
    %% procesing
    % I = imresize(I,[240 320]);
    [M, N] = size(I);
    LPF = avg_filt(avg_filter_size);
    % I1 = I;
    I1 = conv2(I,LPF,'valid');
%     I1 = histeq(I1)
%     I1 = histeq(I1);
    meanGrayLevel = mean2(I1);
    meanGrayLevel = meanGrayLevel/400;
%     if meanGrayLevel > 
    I_bin = ~imbinarize(uint8(I1),meanGrayLevel);
    for i=1:20
        I_bin = medfilt2(I_bin);
    end
    %% matlab nonuniform illumination
    se = strel('disk', STRUCTURING_ELEMENT_SIZE);
    I_bin = imclose(I_bin, se);
    
    
    I_cell = {I, I1, I_bin};
    I_props = regionprops(I_bin);

    digits = cellmat(0);
    j = 1;
    for i=1:length(I_props)
        if I_props(i).Area > (M*N / 25)
%                 rect = rectangle('Position',I_props(i).BoundingBox,...
%                     'EdgeColor','r','LineWidth',3);
                [digit_bin, pre, pos] = I_crop_withBound(I_bin,I_props(i).BoundingBox);
                
                x_origin = I_props(i).BoundingBox(1)-pre(2);
                y_origin = I_props(i).BoundingBox(2)-pre(1);
                
                x_width = I_props(i).BoundingBox(3) + pre(2)+pos(2);
                y_width = I_props(i).BoundingBox(4) + pre(1)+pos(1);
                BoundingBox = [x_origin,y_origin,...
                    x_width, y_width];
                
%                 digitBin = imresize(digit_bin, [digit_side digit_side]);
                digits_bin{j} = digit_bin;

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

% function filt = avg_filt(n)
%     filt = 1/(n^2)*ones(n);
% end
