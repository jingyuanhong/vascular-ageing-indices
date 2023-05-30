load('PW_indx.mat')
clear curr_indx
curr_indx.age = table2array(PW_indx.v(data.plausibility.plausibility_log,1));
curr_indx.v = table2array(PW_indx.v(data.plausibility.plausibility_log,2:end));
curr_indx.unit = PW_indx.unit(2:end)  ;
curr_indx.name = PW_indx.name(2:end)  ;

%% Aortic E
seg_no = [1];%,2,14,18,27];%,28,35,37,39,41];%segment number 
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
E_modulus = sum((rho_weight.*((20/3)*((k(:,1).*exp(k(:,2).*radius_cm))+k(:,3))/10)),2); %g*s-2*cm-1 -> kg*s-2*m-1
%% Bland-Altman Plot
ref = sqrt(E_modulus/(10*1060));
% ref = all_var.v(:,1);
tst = curr_indx.v(:,1:6);
age = struct2table(data.haemods).age(data.plausibility.plausibility_log);
paper_size = [900, 900];
fig_settings.lwidth = 1.2;
fig_settings.ftsize = 10;
fig_settings.color_s = 0*[1 1 1];
fig_settings.color_d = 0.2*[1 1 1];
figure('Position', [20,20,paper_size])
curr_color = {[0 33 245]/255 [139 143 247]/255 [204 205 251]/255 [247 206 205]/255 [240 146 143]/255 [235 51 35]/255};
labels = {'25 yo','35 yo','45 yo','55 yo','65 yo','75 yo'};
tiledlayout(3,2,'TileSpacing','Compact');
CI = zeros(1,width(tst));
ae = zeros(1,width(tst));
sd = zeros(2,width(tst));
for i  = 1:width(tst)
    dif = ref - tst(:,i);
    ave = (ref + tst(:,i))/2;
    
    
    mean_dif = mean(dif);
    std_dif = std(dif);
    mid_line = 2:0.001:18;
    ae(i) = mean_dif;
    sd(1,i) = mean_dif+1.96*std_dif;
    sd(2,i) = mean_dif-1.96*std_dif;
%     fig{i} = subplot(1,3,i);
    fig{i} = nexttile;
    for age_lv = 1:6
        dif_aged = dif(age == ((age_lv*2-1)+4)*5);
        ave_aged = ave(age == ((age_lv*2-1)+4)*5);
        plot(fig{i},ave_aged,dif_aged,'o','Color',curr_color{age_lv},'MarkerSize', 2)
        hold on
    end
    plot(fig{i},mid_line,mean_dif*ones(1,length(mid_line)),'-','Color',fig_settings.color_s, 'LineWidth', fig_settings.lwidth)
    hold on 
    plot(fig{i},mid_line,(mean_dif+1.96*std_dif)*ones(1,length(mid_line)),':','Color',fig_settings.color_d, 'LineWidth', fig_settings.lwidth)
    hold on
    plot(fig{i},mid_line,(mean_dif-1.96*std_dif)*ones(1,length(mid_line)),':','Color',fig_settings.color_d, 'LineWidth', fig_settings.lwidth)
    hold on
%     plot(ave,dif,'o','Color',[0.49,0.18,0.56],'MarkerSize', 2)  
    xlim([5,20])
    ylim([-4.1 2.5])
    ylabel(['aoPWV_t - ' curr_indx.name{i}])
    if ismember(i,[5 6])
        xlabel('mean aoPWV_t' )
    else
        xticklabels({})
    end
    tmp = find(dif > (mean_dif+1.96*std_dif) | dif < (mean_dif-1.96*std_dif));
    ci = length(tmp)/length(dif);
    CI(i) = 1-ci;
end
PATH_SAVE = [PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\figure\fig 7\'];
% saveas(gcf, [PATH_SAVE 'BA plot-Theor aoroot PWV_new.png'])
% clearvars -except data PATH_ROOT
%%
a = find(ave <13.6720 & ave >13.6719);