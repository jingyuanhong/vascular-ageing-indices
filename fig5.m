load('PW_indx.mat')
clear curr_indx
curr_indx.age = table2array(PW_indx.v(data.plausibility.plausibility_log,1));
curr_indx.v = table2array(PW_indx.v(data.plausibility.plausibility_log,2:end));
curr_indx.unit = PW_indx.unit(2:end)  ;
curr_indx.name = PW_indx.name(2:end)  ;
%% Aveaged E
seg_no = [1,2,14,18,27,28,35,37,39,41];%segment number 
%'asc_aorta''aortic_arch1''aortic_arch2''desc_thor_aorta1''desc_thor_aorta2''abd_aorta1''abd_aorta2''abd_aorta3''abd_aorta4''abd_aorta5'
sub = 1:length(data.haemods);
art_len = zeros(length(sub(data.plausibility.plausibility_log)),length(seg_no)); %arterial length
rho_weight = zeros(length(sub(data.plausibility.plausibility_log)),length(seg_no)); %weights for young's modulus calculation
j = 1;
for i = sub(data.plausibility.plausibility_log)%1:length(data.haemods)
    for p = 1:length(seg_no)
        art_pos = find(cell2mat(data.path_waves.aorta_foot(i).segment_no) == seg_no(p));%arterial position
        art_len(j,p) = sum(data.path_waves.aorta_foot(i).artery_dist(art_pos));
    end
rho_weight(j,:) = art_len(j,:)/sum(art_len(j,:));
j = j+1;
end
radius = (data.config.network.inlet_radius(data.plausibility.plausibility_log,seg_no)+data.config.network.outlet_radius(data.plausibility.plausibility_log,seg_no))/2;
radius_cm = radius*100;
k = data.config.constants.k(data.plausibility.plausibility_log,:);

% PATH_SAVE = [PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\figure\figure 3\central\'];
E_modulus = sum((rho_weight.*((20/3)*((k(:,1).*exp(k(:,2).*radius_cm))+k(:,3))/10)),2);
clear seg_no art_len rho_weight radius radius_cm
PATH_SAVE = [PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\figure\fig 5\'];
% fig_plot(f3_var,E_modulus,PATH_SAVE)

rel_rsq_s.all = [];
rel_rsq_p.all = [];

rel_rsq_s.sub25 = [];
rel_rsq_p.sub25 = [];

rel_rsq_s.sub75 = [];
rel_rsq_p.sub75 = [];
%%
[rel_rsq_s,rel_rsq_p] = fig_plot_bi(curr_indx,E_modulus,PATH_SAVE,rel_rsq_s,rel_rsq_p);
% clearvars -except data PATH_ROOT
%%
function [rel_rsq_s,rel_rsq_p] = fig_plot_bi(Var,E,PATH_SAVE,rel_rsq_s,rel_rsq_p)
    if nargin < 3
        PATH_SAVE = [PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\figure\fig 5\'];
    end
    var_arry = Var.v;
    var_name = Var.name;
    var_unit = Var.unit;
    var_age = Var.age;
    
    age_color = {[0 33 245]/255  [235 51 35]/255};
    E_M = E./1e6; % Pa->MPa
    paper_size = [650,800];
    ftsize = 9;
    figure('Position', [20,20,paper_size]);
    % [2 3 4 5 6 7 9 10 11] for PWV ylims -> [4.2,20]
    tiledlayout(5,3,'TileSpacing','Compact');
    for i = 1:width(var_arry)%[3 4 2 5 6 7 9 10 11 8 12 13 14 15 16 17 20 21 22 23 24 25 26 27] %selected parameters
%         no = [3,5,9,10,13,12,14,22,17,24,26,8];
        
        nexttile
        plot(E_M,var_arry(:,i), '.', 'MarkerSize', 6, 'Color', 0.2*[1 1 1]);
        hold on
    %     h = lsline; %least-squares line
    %     h.Color = 'k';
 
        age_25 = find(var_age == 25);
        plot(E_M(age_25),var_arry(age_25,i), '.', 'MarkerSize', 6, 'Color',age_color{1});
        hold on
        age_75 = find(var_age == 75);
        plot(E_M(age_75),var_arry(age_75,i), '.', 'MarkerSize', 6, 'Color',age_color{2});

        
        set(gca, 'FontSize', ftsize,'Box','on','LineWidth', 0.5)
       if ismember(i,[13 14 15])
            xlabel('E_{Ao} [MPa]', 'FontSize', ftsize)
        elseif ismember(i,[1 2 3])
            xlabel('E_{Ao} [MPa]', 'FontSize', ftsize)
            set(gca,'XAxisLocation','top')
        else
            xticklabels({})
        end
%         xlim([0,35])
        if ismember(i,[1 2 3 4 5 6])
            ylim([4.2,20])
        end
        ylabel([var_name{i} var_unit{i}], 'FontSize', ftsize)
    %     xticklabels({})
        refs = E_M;
        vals = var_arry(:,i);
        
        rel_els = find(~isnan(vals));
        temp = corr(refs(rel_els), vals(rel_els),'Type','Spearman');
        rel_rsq_s.all(i) = temp;
        clear temp rel_els
    
        rel_els = find(~isnan(vals));
        temp = corr(refs(rel_els), vals(rel_els),'Type','Pearson');
        rel_rsq_p.all(i) = temp;
        clear temp rel_els
        
        refs = E_M(age_25);
        vals = var_arry(age_25,i);
        
        rel_els = find(~isnan(vals));
        temp = corr(refs(rel_els), vals(rel_els),'Type','Spearman');
        rel_rsq_s.sub25(i) = temp;
        clear temp rel_els
    
        rel_els = find(~isnan(vals));
        temp = corr(refs(rel_els), vals(rel_els),'Type','Pearson');
        rel_rsq_p.sub25(i) = temp;
        clear temp rel_els
        
        refs = E_M(age_75);
        vals = var_arry(age_75,i);
        
        rel_els = find(~isnan(vals));
        temp = corr(refs(rel_els), vals(rel_els),'Type','Spearman');
        rel_rsq_s.sub75(i) = temp;
        clear temp rel_els
    
        rel_els = find(~isnan(vals));
        temp = corr(refs(rel_els), vals(rel_els),'Type','Pearson');
        rel_rsq_p.sub75(i) = temp;
        clear temp rel_els
        
    end
    rel_rsq_s.all = rel_rsq_s.all';
    rel_rsq_p.all = rel_rsq_p.all';
    rel_rsq_s.sub25 = rel_rsq_s.sub25';
    rel_rsq_p.sub25 = rel_rsq_p.sub25';
    rel_rsq_s.sub75 = rel_rsq_s.sub75';
    rel_rsq_p.sub75 = rel_rsq_p.sub75';
%     saveas(gcf, [PATH_SAVE 'E_vs_indx.png'])
%         close all;
end