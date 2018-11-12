clc, clear, close all
addpath('VLFEATROOT')

run('VLFEATROOT/toolbox/vl_setup');
%vl_version verbose
%for the closing operation
AVG_FILTER_SIZE = 3; %for extracting digits
DIGIT_SIZE = 120; %image size
addpath('numbers'); %where images/mat files are stored are stored
%% Loading and pre-processing images
%image A & B
load('digits_map.mat');
% scores_collected = [];
% no_matches = [];
for i=0:9
    no_matches = [];
    scores_collected = [];

    dig = digits_map(i);
    Ia_orig = dig{3};    
    % extract 120 x 120 digits normal and binarized
    [Ia, Ia_bin] = extract_digits(Ia_orig,AVG_FILTER_SIZE,DIGIT_SIZE);
    Ia_bin = Ia_bin{1};
    Ia = process_I(Ia{1});
    dig2 = digits_map(i);
    Ib_orig = dig2{5};
%     Ia = Ia{1};
    [Ib, Ib_bin] = extract_digits(Ib_orig, AVG_FILTER_SIZE, DIGIT_SIZE);
%         Ib = process_I(Ib{1});
%     Ib = Ib{1};
    % pre processing (adding blur to binarized and getting canny edges
    Ia_bin = process_bin_num(Ia_bin,2);
    Ib_bin = process_bin_num(Ib_bin{1},2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SIFT

%         [fa,da] = run_SIFT(Ia,0);
%         [fb,db] = run_SIFT(Ib, size(Ia,2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [fa_bin,da_bin] = run_SIFT(Ia_bin,0);
    x_offset = size(Ia_bin,2);
    [fb_bin,db_bin] = run_SIFT(Ib_bin,x_offset);
    %%  SIFT MATCHES ///////////////////////////////////////////////////////////
    % --------------------------- grayscale
%         [matches_normal, scores_normal] = vl_ubcmatch(da, db,1) ;
%         m  = median(scores_normal);
%         std_n = std(scores_normal);
%         s = 0; % .5*std_n;
%         index_m = find(scores_normal < (m +std_n));
%         selec_scores_norm = scores_normal(index_m);
%         selec_matches_norm = matches_normal(:,index_m);

    % --------------------------- binarized
    [matches_bin, scores_bin] = vl_ubcmatch(da_bin, db_bin,1.3);
    
    [selec_matches, selec_scores ] = nonUniques(matches_bin, scores_bin);
    
%     m_bin  = median(scores_bin);
%     std_bin = std(scores_bin);
%     index_m_bin = find(scores_bin < m_bin +1.5*std_bin);
%     selec_scores = scores_bin(index_m_bin);
%     selec_matches = matches_bin(:,index_m_bin);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         figure
%         subplot(121);
%         imshow(uint8([Ia,Ib]));
%         title('Ia X  Ib');
%         draw_lines(da,fa,db,fb,selec_matches_norm,selec_scores_norm,x_offset);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        subplot(2,5,i+1);
%         imshow(uint8([Ia,Ib{1}]));
        imshow([Ia_bin,Ib_bin]);
        title('Ia bin X Ib bin');
% 
        draw_lines(da_bin, fa_bin, db_bin, fb_bin,...
            selec_matches, selec_scores, x_offset);
end
%%

function [unique_maches, unique_scores] = nonUniques(matchings,scorings)
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

function draw_lines(da,fa,db,fb,matches,scores,x_offset)
    xa = fa(1,matches(1,:)) ;
    xb = fb(1,matches(2,:)) + x_offset;
    ya = fa(2,matches(1,:)) ;
    yb = fb(2,matches(2,:)) ;
     
    hold on ;
    h = line([xa ; xb], [ya ; yb]) ;
    set(h,'linewidth', 1, 'color', 'b') ;

    vl_plotframe(fa(:,matches(1,:)));
    fb(1,:) = fb(1,:) + x_offset;
    vl_plotframe(fb(:,matches(2,:)));
    colormap gray
    axis image off ;
    hold off

end

function I = process_bin_num(I_bin,avg_filter_size)
    LPF = avg_filt(avg_filter_size);
    I_bin = uint8(I_bin*255);
   I = I_bin;
    % I1 = I;
%     I = conv2(I_bin,LPF,'valid');
    I = edge(I,'Canny');
    x = fspecial('motion',10,0);
%     y = fspecial('motion',25,90);
%     LPF = avg_filt(5);
%     
     I = conv2(I,x,'valid');
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