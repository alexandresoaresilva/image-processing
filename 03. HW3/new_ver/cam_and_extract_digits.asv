clc, clear, close all
addpath('VLFEATROOT')

run('VLFEATROOT/toolbox/vl_setup');
addpath('numbers'); %where images/mat files are stored are stored



%% paramteters
AVG_FILTER_SIZE = 3; %for extracting digits
DIGIT_SIZE = 120; %Size of square-cropped digit
MOTION = 10;      %Experimenting with motion blur (adding predictable complexity to character)
MATCH_THRESHOLD = 1.5;  %
BINARIZE = 0;
SHOW_MATCHES = 1;

%% Loading and pre-processing images
%image A & B
load('digits_map.mat');
dig = digits_map(5);
Ia_orig = dig{3};
% [Ia, Ia_bin] = extract_digits(Ia_orig,AVG_FILTER_SIZE,DIGIT_SIZE);
% 
% if BINARIZE
%     Ia = process_bin_num(Ia_bin{1}, 3, MOTION );
% else
%     Ia = Ia{1};
% end

% % figure; clf;
% for i=0:9
%     no_matches = [];
%     scores_collected = [];
%     % extract 120 x 120 digits normal and binarized
%     %digit to compare
%     dig2 = digits_map(i);
%     Ib_orig = dig2{5};
%     [Ib, Ib_bin] = extract_digits(Ib_orig, AVG_FILTER_SIZE, DIGIT_SIZE);
%     
%     if BINARIZE
%         Ib = process_bin_num(Ib_bin{1}, 3, MOTION );
%     else
%         Ib = Ib{1};
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% SIFT
%     
%     [fa_bin,da_bin] = vl_sift(single(Ia)) ;
%     [fb_bin,db_bin] = vl_sift(single(Ib)) ;
%     %%  SIFT MATCHES ///////////////////////////////////////////////////////////
%     % --------------------------- binarized
%     [matches_bin, scores_bin] = vl_ubcmatch(da_bin, db_bin, MATCH_THRESHOLD);
%     [matches_bin, scores_bin] = remove_inf_and_nans(matches_bin, scores_bin);
%     [matches_bin,scores_bin] = unique_matches_scores(fa_bin, 1, matches_bin,scores_bin);
%     [selec_matches,selec_scores] = unique_matches_scores(fb_bin, 0, matches_bin, scores_bin);
%     [selec_matches, selec_scores] = remove_outliers(selec_matches, selec_scores, 0);
%     [selec_matches, selec_scores]  = remove_senseless_scores(selec_matches, selec_scores, fa_bin, fb_bin);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     subplot(2,5,i+1);
% %     imshow([Ia_bin,Ib_bin]);
% %     title('Ia bin X Ib bin');
% %     plot_SIFT_lines(Ia, Ib, fa_bin, fb_bin, selec_matches);
% end
% [most_freq_num, scores_sum,no_matches, total_matches] = match_and_predict(Ia_orig, digits_map, SHOW_MATCHES,...
%     AVG_FILTER_SIZE, MATCH_THRESHOLD, DIGIT_SIZE, BINARIZE ,MOTION)

[predictedNumberOptions, minScores, scores, normalized_scores,...
    no_matches, counts] =...
    matchAndPredict(digits_map, Ia_orig, 0, 1.5, AVG_FILTER_SIZE, DIGIT_SIZE, 1);


function [most_freq_num, scores_sum,no_matches, total_matches] =...
    match_and_predict(I, digits_map,show_matches, avg_filter_size,...
                      match_thresh, digit_size, binarize, motion_size)
    
    [Ia, Ia_bin] = extract_digits(I,avg_filter_size,digit_size);
    if binarize
       Ia = process_bin_num(Ia_bin{1},2, motion_size);
    else
       Ia = Ia{1};
    end
    [fa,da] = vl_sift(single(Ia)) ;
    scores_sum  = []; 
    no_matches = [];
    total_matches = [];
    %looks for the digit match
    for i=1:length(digits_map)
        test_dig = digits_map(i-1);
        figure;
        for j=1:length(test_dig)
           [Ib, Ib_bin] =...
               extract_digits(test_dig{j}, avg_filter_size, digit_size);
           if binarize %process it in binary with motion blur
               Ib = process_bin_num(Ib_bin{1},2, motion_size);
           else
               Ib = Ib{1};
           end
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %% SIFT
            [fb,db] = vl_sift(single(Ib));
            [matches, scores] = vl_ubcmatch(da, db, match_thresh);
            [matches, scores] = remove_inf_and_nans(matches, scores);
            [matches,scores] = unique_matches_scores(fb,0, matches, scores);    
            [matches,scores] = unique_matches_scores(fa,1,matches, scores);    
            [matches, scores] = remove_senseless_scores(matches, scores, fa, fb);
            [matches,scores] = keep_only_three_best_features(matches, scores);
            
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
    [~,i] = min(sum(scores_sum,2));
%     most_freq_num = i-1;
     most_freq_num = find_match(scores_sum, total_matches);
%     most_freq_num = find_most_frequent_digit(scores_sum);
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

function [matches,scores] = unique_matches_scores(fb_or_fa, is_fa, matches,scores)
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
function [predictedNumberOptions, minScores, scores,...
            normalized_scores, no_matches, counts] = ...
        matchAndPredict(templateImageMap, Ia,... %function definition
            showMatches, SIFT_match_thresh, avg_filter_size,...
            digit_size, binarized)
        
    digit = templateImageMap;
    [fa,da] = vl_sift(im2single((Ia))) ;
    
    for j=1:length(digit)
        digit_no = digit(j-1);
        
        for i=1:length(digit_no)
            Ib_orig = digit_no{i};
%             [~, Ibin] = extract_digits(Ib_orig, avg_filter_size, digit_size);
%             Ib = process_bin_num(Ib{1},2);
            [Ib, Ibin] = extract_digits(Ib_orig, avg_filter_size, digit_size);
            
            if binarized
              	Ib = process_bin_num(Ibin{1},2, 15);
            else
                Ib = Ib{1};
            end
            
            [fb,db] = vl_sift(im2single((Ib))) ;
            [matches, score] = vl_ubcmatch(da,db,SIFT_match_thresh);
%             [matches, score] = clean_inf_and_nans(matches, score);
            [matches, score] = nonUniqueReduction(matches,score);
%             [matches,score] = unique_matches_scores(fb, matches, score);
%             [matches, score]  = remove_senseless_scores(matches, score, fa, fb);
            if(showMatches)
                figure; clf ;
                plot_SIFT_lines(Ia, Ib, fa, fb, matches);
            end

            if(length(score)>4)
                [matches, score] =...
                    ignore_scores_above_median(matches, score, 0);
            end
            %colecting the sum of scores to find the smallest among the
            %most frequent 
            scores(j,i) = sum(score);
            no_matches(j,i) = length(matches);
        end
    end
    
    for k=1:size(scores,2)
        normalized_scores(:,k) = abs((scores(:,k)-mean(scores(:,k)))./no_matches(:,k));
    end
    
    [minScores, predictedNumberOptions] = min(normalized_scores); 
    nans = isnan(minScores);
    minScores(nans) = inf;
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


function [matches, scores] = remove_senseless_scores(matches, scores, fa, fb)


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
        xb = fb(1,matches(2,:)) ;
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
    
    [matches, scores] = remove_outliers(matches, scores, 0);
end

function [unique_maches, unique_scores] = nonUniqueReduction(matchings,scorings)
    Aa = matchings(1,:);
    Ab = matchings(2,:);
    Bc = scorings;
    if unique(scorings) == 0 %Don't run function if input image is identical to a sample image
        unique_maches = matchings;  %i.e. every match is already perfect
        unique_scores = scorings;   
    else
        n = length(Bc);         %Sort scores from lowest to highest
        while n ~= 0            %and align corresponding indices for matches
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
        for j = 1:length(B2)           %Prioritizing closest match, reduce non-uniqe matches
            P = A2(1,j);               %Feature index of first image
            if(isempty(find(unique_maches == P,1)))
                Q = A2(2,j);           %Feature index of second image
                idxOfUniqueMinScore = find(B2 == min(B2(round((find(A2 == Q).')/2)))); %Find all corresponding matches
                                                                                        %and select match with minimum score
                if(isempty(find(unique_maches == Q,1)))
                    unique_maches(1,k) = A2(1,j);                       %If match pair doesn't already exist
                    unique_maches(2,k) = A2(2,idxOfUniqueMinScore);     %Add to reduced match and score set
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
    disp(['Real: ',num2str(i), ' Predicted Number: ', num2str(pI-1)])   %Display two best predictions
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
    LPF = avg_filt(avg_filter_size);
    I_bin = uint8(I_bin*255);
   I = I_bin;
    % I1 = I;
    I = conv2(I_bin,LPF,'valid');
%     I = edge(I,'Canny');
%     x = fspecial('motion',motion_blur_amount,0);
%     y = fspecial('motion',25,90);
%     LPF = avg_filt(5);
%     
%      I = conv2(I,x,'valid');
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