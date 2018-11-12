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
    Ia_orig = dig{1};    
    % extract 120 x 120 digits normal and binarized
    [~, Ia_bin] = extract_digits(Ia_orig,AVG_FILTER_SIZE,DIGIT_SIZE);
    Ia_bin = Ia_bin{1};
    x_offset = size(Ia_bin,2);

    for j=0:9
        dig2 = digits_map(j);
        Ib_orig = dig2{1};
    %     Ia = Ia{1};
        [Ib, Ib_bin] = extract_digits(Ib_orig, AVG_FILTER_SIZE, DIGIT_SIZE);
%         Ib = process_I(Ib{1});
    %     Ib = Ib{1};
        % pre processing (adding blur to binarized and getting canny edges
        Ia_bin = process_bin_num(Ia_bin,2);
        Ib_bin = process_bin_num(Ib_bin{1},2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% SIFT
        % --------------------------- binarized
        [fa_bin,da_bin] = run_SIFT(Ia_bin,0);
        [fb_bin,db_bin] = run_SIFT(Ib_bin,x_offset);
        %%  SIFT MATCHES ///////////////////////////////////////////////////////////
        % --------------------------- binarized
        [matches_bin, scores_bin] = vl_ubcmatch(da_bin, db_bin) ;
        
        m_bin  = median(scores_bin);
        std_bin = std(scores_bin);
        index_m_bin = find(scores_bin < m_bin +1.5*std_bin);
        selec_scores = scores_bin(index_m_bin);
        selec_matches = matches_bin(:,index_m_bin);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%       figure
%         imshow([Ia_bin,Ib_bin]);
%         title('Ia bin X Ib bin');
% 
%         draw_lines(da_bin, fa_bin, db_bin, fb_bin,...
%             selec_matches, selec_scores, x_offset);

        no_matches = length(selec_matches);
        scores_collected(j+1) = sum(selec_scores)./no_matches;
    end
    clear Ia;
    predict_no(i, selec_matches, selec_scores);
end
%% 
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
   
    % I1 = I;
    I = conv2(I_bin,LPF,'same');
    I = edge(I,'Canny');
    x = fspecial('motion',5,0);
%     y = fspecial('motion',25,90);
%     LPF = avg_filt(5);
%     
     I = conv2(I,x,'same');
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