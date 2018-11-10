function sub_plots(I_cell)
    n = length(I_cell);
    no_of_plots_per_row = ceil(n/2);
    figure;
    for i=1:n
        subplot(2,no_of_plots_per_row,i)
        imshow(I_cell{i},[]);
        title(num2str(i));
    end
end