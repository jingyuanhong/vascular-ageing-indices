load('PW_indx.mat')
clear curr_indx
curr_indx.age = table2array(PW_indx.v(:,1));
curr_indx.v = table2array(PW_indx.v(:,2:end));
curr_indx.unit = PW_indx.unit(2:end);
curr_indx.name = PW_indx.name(2:end);

%%
all_params = data.config.variations.params;
param_names = data.config.variations.param_names;
[rel_sims,  baseline_logs, baseline_sims, param_variations, param_values] = deal(cell(length(param_names),1));

%%
%cycle through each simulation input parameter
paper_size = [1000, 1000];
figure('Position', [20,20,paper_size])
tiledlayout(5,3,'TileSpacing','Compact')
i_vals = zeros(width(curr_indx.v),19);
for i = 1:width(curr_indx.v)
for param_no =  1 : length(param_names)
    
    % identify those columns which correspond to other parameters (i.e. not this one)
    columns_for_other_params = setxor(1:length(param_names), param_no); 
    % extract variations corresponding to the other parameters
    temp = all_params(:,columns_for_other_params);
    % identify simulations in which none of these other parameters were varied from baseline values 
    rel_sims{param_no} = find(~any(temp,2));
    % skip this parameter if it was not varied at all from its baseline value
    if length(rel_sims{param_no}) <= sum(data.config.baseline_sim_for_age)
        continue
    end % select independent index
    
    % use all the simulations
    rel_sims{param_no} = 1:length(data.config.age); % subject number
    % setup variables
    temp_baseline_vals = nan(length(rel_sims{param_no}),1);
    temp_baseline_logs = false(length(rel_sims{param_no}),1);
    % cycle through each simulation
    for s = 1 : length(rel_sims{param_no})
        % identify age of this simulation
        curr_sim = rel_sims{param_no}(s);
        curr_age = data.config.age(curr_sim);
        % extract baseline value of this parameter
        temp_rel_baseline_age_sim = find(data.config.age == curr_age & data.config.baseline_sim_for_age);
        eval(['temp_baseline_vals(s) = data.config.' param_names{param_no} '(temp_rel_baseline_age_sim);']);
        % store the baseline simulation for this age
        if temp_rel_baseline_age_sim == curr_sim
            temp_baseline_logs(s) = true;
        end
        temp_baseline_sims(s) = temp_rel_baseline_age_sim;
    end
    % store details of the baseline simulations corresponding to each simulation
    baseline_logs{param_no} = temp_baseline_logs; clear temp_baseline_logs
    baseline_sims{param_no} = temp_baseline_sims; clear temp_baseline_sims
    baseline_vals{param_no} = temp_baseline_vals; clear temp_baseline_vals
    % store the variations for this parameter (in number of SDs from age-specific mean)
    param_variations{param_no} = all_params(rel_sims{param_no},param_no);
    % store the values of this parameter
    eval(['param_values{param_no} = data.config.' param_names{param_no} '(rel_sims{param_no});']);
end
%% sensitivity analyses

    I = nan(length(param_names), 1);
    for param_no = 1 : length(param_names)
        % skip if this parameter wasn't varied indpendently
        if length(param_variations{param_no}) <= 1 % (1 accounts for the baseline sim)
            I(param_no,:) = nan;
            continue
        end
        
        % extract data for this parameter
%         curr_rel_sims = rel_sims{param_no}(data.plausibility.plausibility_log);
%         curr_rel_variations = param_variations{param_no}(data.plausibility.plausibility_log);
%         curr_baseline_vals = baseline_vals{param_no};
%         curr_baseline_sim_nos = baseline_sims{param_no};
%         curr_baseline_logs = baseline_logs{param_no};
%         
%         % perform sensitivity analysis
%         indx_v = curr_indx.v(:,i);
% 
%         % identify relevant simulations (excluding those which were physiologically implausible)
%         % rel_els = ~curr_baseline_logs & curr_rel_variations~=0 & ~isnan(all_values) & data.plausibility.plausibility_log;
%         rel_els = ~curr_baseline_logs(data.plausibility.plausibility_log) & curr_rel_variations~=0 & ~isnan(indx_v(data.plausibility.plausibility_log)) & data.plausibility.plausibility_log(data.plausibility.plausibility_log);
% 
%         rel_value_sims = curr_rel_sims(rel_els);
%         rel_baseline_sims = curr_baseline_sim_nos(rel_els); % each age level
%         rel_variations = curr_rel_variations(rel_els);
%         
%         % identify param values
%         res_param_values = indx_v(rel_value_sims);
%         res_param_baseline_values = indx_v(rel_baseline_sims);
%         
%         % calculate index
%         v = rel_variations;
%         I(param_no) = 100*mean((res_param_values - res_param_baseline_values)./(abs(res_param_baseline_values).*v));
        curr_rel_sims = rel_sims{param_no};
        curr_rel_variations = param_variations{param_no};
        curr_baseline_vals = baseline_vals{param_no};
        curr_baseline_sim_nos = baseline_sims{param_no};
        curr_baseline_logs = baseline_logs{param_no};
        
        % perform sensitivity analysis
        indx_v = curr_indx.v(:,i);

        % identify relevant simulations (excluding those which were physiologically implausible)
        rel_els = ~curr_baseline_logs & curr_rel_variations~=0 & ~isnan(indx_v) & data.plausibility.plausibility_log;

        rel_value_sims = curr_rel_sims(rel_els);
        rel_baseline_sims = curr_baseline_sim_nos(rel_els); % each age level
        rel_variations = curr_rel_variations(rel_els);
        % identify param values
        res_param_values = indx_v(rel_value_sims);
        res_param_baseline_values = indx_v(rel_baseline_sims);
        v = rel_variations;

%         % identify relevent simulations of single input parameter variation
%         % with other input parameters keeping in constants
%         
%         rel_ind_value_sims = rel_value_sims(rel_value_sims<=78);
%         rel_ind_baseline_sims = rel_baseline_sims(rel_value_sims<=78);
%         rel_ind_variations = rel_variations(rel_value_sims<=78);
%         
%         % identify param values
%         res_param_values = indx_v(rel_ind_value_sims);
%         res_param_baseline_values = indx_v(rel_ind_baseline_sims);
%         v = rel_ind_variations;
        % calculate index
        
        I(param_no) = 100*mean((res_param_values - res_param_baseline_values)./(abs(res_param_baseline_values).*v));

    end
%% Make sensitivity plots for each resultant parameter
%     paper_size = [200,200,600,400];

    
    curr_i_vals = I;
    curr_param_names = param_names;
    rel_els = ~strcmp(curr_param_names, 'ht') & ~strcmp(curr_param_names, 'rho');
    curr_i_vals = curr_i_vals(rel_els);
    curr_param_names = curr_param_names(rel_els);
    %req_order = {'hr','sv','lvet','t_pf','reg_vol','mbp','pwv','p_out','dia','len','pvc'};
    
    ref_param = curr_indx.name{i};%resultant_params{res_param_no};
    curr_param_names = strrep(curr_param_names, 'reg_vol', 'Rvol');
    curr_param_names = strrep(curr_param_names, 'p_out', 'Pout');
    
    % ignore parameters which we're not interested in
    rel_no = strcmp(curr_param_names, 'dbp');
    curr_i_vals(rel_no) = nan;
    i_vals(i,:) = curr_i_vals';
    nexttile
    fig_done = make_sens_analysis_fig(curr_param_names, curr_i_vals, 'norm', paper_size, ref_param,i);

end
PATH_SAVE = [PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\figure\fig 6\'];
% saveas(gcf, [PATH_SAVE 'sensIndex_all.png'])
% clearvars -except data PATH_ROOT
%%
function fig_done = make_sens_analysis_fig(model_params, i_vals, type, paper_size, ref_param,sim_no)

% exclude variables which weren't varied
rel_els = ~isnan(i_vals);
i_vals = i_vals(rel_els);
model_params = model_params(rel_els);

if isempty(model_params)
    fig_done = 0;
    return
else
    fig_done = 1;
end

% re-arrange

rel_order = {'hr', 'sv', 'lvet', 'dia', 'pwv', 'mbp'};

for s = 1 : length(rel_order)
    order(s) = find(strcmp(model_params,rel_order{s}));    
end
model_params = model_params(order);
i_vals = i_vals(order);

% re-name
model_params = strrep(model_params, 'hr', 'Heart Rate');
model_params = strrep(model_params, 'sv', 'Stroke Volume');
model_params = strrep(model_params, 'lvet', 'Duration Systole');
model_params = strrep(model_params, 't_pf', 'Time to Peak Flow');
model_params = strrep(model_params, 'Rvol', 'Regurgitation Volume');
model_params = strrep(model_params, 'dia', 'Large Art. Diameter');
model_params = strrep(model_params, 'len', 'Prox. Aortic Length');
model_params = strrep(model_params, 'pwv', 'Input PWV');
model_params = strrep(model_params, 'mbp', 'Periph. Vasc. Res.');
model_params = strrep(model_params, 'pvc', 'Periph. Vasc. Comp.');
model_params = strrep(model_params, 'Pout', 'Outflow Pressure');
model_params = strrep(model_params, 'pvr', 'Periph. Vasc. Res.');   

if strcmp(type, 'abs')
    ylims = [-30, 75];
    i_vals = abs(i_vals);
    ylab_txt = 'abs(I) [%]';
else
    ylims = [-50, 67];
    ylab_txt = 'I [%]';
end

model_params = strrep(model_params, '_', ' ');

ftsize = 10;
fig = bar(i_vals,'EdgeColor',0.4*[1 1 1],'FaceColor',0.4*[1 1 1]);
if ismember(sim_no,[13 14 15])
    set(gca,'XTickLabel',model_params)
else
    set(gca,'XTickLabel',[])
end
set(gca, 'FontSize', ftsize)
if ismember(sim_no,[1 4 7 10 13])
    ylabel(ylab_txt, 'FontSize', ftsize);
elseif ismember(sim_no,[3 6 9 12 15])
    ylabel(ylab_txt, 'FontSize', ftsize);
    set(gca,'YAxisLocation','right')
else
    set(gca,'YTickLabel',[])
end

% ylab = ylabel(ylab_txt, 'Rotation', 0, 'FontSize', ftsize);
% set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
%ylim(ylims)
% [label, units, abbr, graph_title, graph_title_no_units] = make_param_label(ref_param);
title(ref_param, 'FontSize', ftsize)
xtickangle(30);
xlim([0.5, length(model_params)+0.5])

hold on
% ylims = [min(i_vals) - 0.1*range(i_vals), max(i_vals) + 0.22*range(i_vals)];
% ylims = [-15 20];
plot(3.5*ones(1,2), [-100 100], '--k')
ylim(ylims)

end