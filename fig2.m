clear all;clc
% PATH_ROOT = 'C:\Users\43897\';
PATH_ROOT = 'C:\Users\Jingyuan Hong\';
% PATH_ROOT = 'C:\Users\jh22\';
load([PATH_ROOT,'Downloads\pwdb_data_w_aorta_finger_path'])
%%
clear rel_sims
rel_sims.ages = find(data.config.baseline_sim_for_age);
baseline_age = data.config.age(data.config.baseline_sim_for_all);

% Setup figure for wave speed curves
paper_size = [600, 230];
fig_settings.lwidth = 1;
fig_settings.ftsize = 10;
% fig_settings.colors = [1,1,1];

% cycle through each plot type
plot_types = fieldnames(rel_sims);
figure('Position', [20,20,paper_size])
labels = {'25 y/o','35 y/o','45 y/o','55 y/o','65 y/o','75 y/o'};
% cycle through each relevant simulation
eval(['rel_sims.curr = rel_sims.' plot_types{1} ';'])
age = struct2table(data.haemods).age(data.plausibility.plausibility_log,:);
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
E_ao = sum((rho_weight.*((20*((2*k(:,1).*exp(k(:,2).*radius_cm))+2*k(:,3))/3)/1e7)),2);
E_ao_age = [];
for sim_no = 1 : length(rel_sims.curr)
    curr_sim_no = rel_sims.curr(sim_no);

    rel_network_spec.inlet_radius = data.config.network.inlet_radius(curr_sim_no,:);
    rel_network_spec.outlet_radius = data.config.network.outlet_radius(curr_sim_no,:);
    k = data.config.constants.k(curr_sim_no,:);
    ave_radius = mean([rel_network_spec.inlet_radius(:), rel_network_spec.outlet_radius(:)], 2);
    ave_radius_cm = ave_radius*100;



    E = (20*((2*k(1)*exp(k(2).*ave_radius_cm))+2*k(3))/3)/1e7; % g*s-2*cm-1 -> MPa

    clear k rel_network_spec
    [ave_radius_cm, order] = sort(ave_radius_cm);
    E = E(order);
    curr_color = {[0 33 245]/255 [139 143 247]/255 [204 205 251]/255 [247 206 205]/255 [240 146 143]/255 [235 51 35]/255};%color_intensity*fig_settings.colors;

    f1 = subplot(1,2,1);
    plot(f1,ave_radius_cm, E, 'o-', 'Color', curr_color{sim_no}, 'LineWidth', fig_settings.lwidth, ...
        'MarkerSize', 3, 'MarkerFaceColor', curr_color{sim_no}, 'MarkerEdgeColor', curr_color{sim_no})    
    hold on
    clear wave_speed ave_radius curr_color color_intensity

    yo = ((sim_no*2-1)+4)*5;
    E_ao_age = [E_ao_age; yo, mean( E_ao(age == yo,:)) std(E_ao(age == yo,:))];  
end

f2 = subplot(1,2,2);
plot(f2,E_ao_age(:,1),E_ao_age(:,2),'o-','MarkerSize', 4,'Color', 0*[1 1 1])
hold on
plot(f2,E_ao_age(:,1),E_ao_age(:,2)+E_ao_age(:,3),'.--','MarkerSize', 10,'Color', 0.2*[1 1 1])
hold on
plot(f2,E_ao_age(:,1),E_ao_age(:,2)-E_ao_age(:,3),'.--','MarkerSize', 10,'Color', 0.2*[1 1 1])

ylabel(f1,'E [MPa]', 'FontSize', fig_settings.ftsize);
xlabel(f1,'r [cm]', 'FontSize', fig_settings.ftsize);
ylabel(f2,'Aortic E [MPa]', 'FontSize', fig_settings.ftsize);
xlabel(f2,'Age [years]', 'FontSize', fig_settings.ftsize);
l1 = legend(f1,labels, 'FontSize', fig_settings.ftsize-2); clear labels
l1.Position(2) = 0.53;
%     ylim([f1,f2,f3],[0.02 0.22])
ylim([f1,f2],[0.4 5])
% xlim(f2,[20 80])
PATH_SAVE = [PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\figure\fig 2\'];
% saveas(gcf, [PATH_SAVE 'E_ao_age_new.png'])

