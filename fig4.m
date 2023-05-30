%%
load('PW_indx.mat')
clear insilico_data
no = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
age = table2array(PW_indx.v(data.plausibility.plausibility_log,1));
curr_indx = table2array(PW_indx.v(data.plausibility.plausibility_log,no));
for i = 1:width(curr_indx)
    insilico_data.v{i} = [];
    for yo = 25:10:75
        insilico_data.v{i} = [insilico_data.v{i}; yo, mean(curr_indx(age == yo,i)) std(curr_indx(age == yo,i))];   
    end
end
insilico_data.unit = PW_indx.unit(no)  ;
insilico_data.name = PW_indx.name(no)  ;
%%
TABLE_PATH = [PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\figure\fig 4\in vivo data.xlsx'];
clear invivo_data
invivo_data.v = {};
for i = 1:width(curr_indx)
    res_t = readtable(TABLE_PATH,'Sheet',insilico_data.name{i});
    L = height(res_t);
    invivo_data.v{i} = [table2array(res_t(1:L-2,1)),table2array(res_t(1:L-2,4)),table2array(res_t(1:L-2,5))];
end
invivo_data.unit = insilico_data.unit;
invivo_data.name = insilico_data.name;
%%
PATH_SAVE = [PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\figure\fig 4\'];
paper_size = [900, 950];
ftsize = 9;
figure('Position', [20,20,paper_size])
tiledlayout(8,4,'TileSpacing','Compact')
for i = 1:15 %2:25%length(var_name)
    age_vivo = invivo_data.v{i}(:,1);
    y_vivo = invivo_data.v{i}(:,2);
    err_vivo = invivo_data.v{i}(:,3);

    age_silc = insilico_data.v{i}(:,1);
    y_silc = insilico_data.v{i}(:,2);
    err_silc = insilico_data.v{i}(:,3);
    
    
    nexttile
    plot(age_vivo,y_vivo,'-','MarkerSize', 4,'Color', 0*[1 1 1])
    hold on
    plot(age_vivo,y_vivo+err_vivo,'--','MarkerSize', 10,'Color', 0.2*[1 1 1])
    hold on
    plot(age_vivo,y_vivo-err_vivo,'--','MarkerSize', 10,'Color', 0.2*[1 1 1])

    set(gca, 'FontSize', ftsize,'Box','on','LineWidth', 1)
    
    ylabel([invivo_data.name{i} invivo_data.unit{i}], 'FontSize', ftsize)
    if ismember(i,[ 15])
        xlabel('Age [years]', 'FontSize', ftsize)
    elseif ismember(i,[1 2])
        xlabel('Age [years]', 'FontSize', ftsize)
        set(gca,'XAxisLocation','top')
    else
        xticklabels([])
    end
    xticks([20 40 60 80]);
    xlim([20,80])
    grid on
%     while i<6 
%         ylim([2,30])
%         break
%     end
%     while i >5        
%         ylim([min(min(y_vivo)- 1.5*max(err_vivo),min(y_silc) - 1.5*max(err_silc)) , max(max(y_vivo) + 1.5*max(err_vivo),max(y_silc) + 1.5*max(err_silc))])
%         break
%     end
ylim([min(min(y_vivo)- 1.5*max(err_vivo),min(y_silc) - 1.5*max(err_silc)) , max(max(y_vivo) + 1.5*max(err_vivo),max(y_silc) + 1.5*max(err_silc))])
%     get(h(i),'Position');
%     set(h(i),'Position',[0.1+(mod(i-1,4))*0.21 0.8-floor((i-1)/4)*0.14 0.15 0.12]);
    nexttile
    plot(age_silc,y_silc,'o-','MarkerSize', 4,'Color', 0*[1 1 1])
    hold on
    plot(age_silc,y_silc+err_silc,'.--','MarkerSize', 10,'Color', 0.2*[1 1 1])
    hold on
    plot(age_silc,y_silc-err_silc,'.--','MarkerSize', 10,'Color', 0.2*[1 1 1])

    set(gca, 'FontSize', ftsize,'Box','on','LineWidth', 1)

    
    ylabel([insilico_data.name{i} insilico_data.unit{i}], 'FontSize', ftsize)
    if ismember(i,[ 15])
        xlabel('Age [years]', 'FontSize', ftsize)
    elseif ismember(i,[1 2])
        xlabel('Age [years]', 'FontSize', ftsize)
        set(gca,'XAxisLocation','top')
    else
        xticklabels([])
    end
    xticks([20 40 60 80]);
    xlim([20,80])
    grid on
%     while i<6 
%         ylim([2,30])
%         break
%     end
%     while i >5        
%         ylim([min(min(y_vivo)- 1.5*max(err_vivo),min(y_silc) - 1.5*max(err_silc)) , max(max(y_vivo) + 1.5*max(err_vivo),max(y_silc) + 1.5*max(err_silc))])
%         break
%     end
ylim([min(min(y_vivo)- 1.5*max(err_vivo),min(y_silc) - 1.5*max(err_silc)) , max(max(y_vivo) + 1.5*max(err_vivo),max(y_silc) + 1.5*max(err_silc))])
end
% saveas(gcf, [PATH_SAVE 'invivo-insilico.png'])
% clearvars -except data PATH_ROOT