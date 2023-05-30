
PATH_ROOT = 'C:\Users\Jingyuan Hong\';
PATHS.input_folder = [PATH_ROOT,'/OneDrive - King''s College London/PhD Work/code/PW_studies/'];
PATHS.output_folder = [PATH_ROOT,'/OneDrive - King''s College London/PhD Work/2022/Comprehensive Vascular Ageing Analysis/figure/'];
PATHS.exported_data = [PATHS.input_folder, 'pwdb_' num2str(pwdb_no), filesep, 'exported_data', filesep];
PATHS.Analysis_figures = [PATHS.output_folder, 'figure 1/'];
PATHS.baseline_waves_fig_a = [PATHS.Analysis_figures, 'baseline_waves_fig_a'];
PATHS.baseline_waves_fig_b = [PATHS.Analysis_figures, 'baseline_waves_fig_b'];
load(PATHS.exported_data_mat_pwdb_data)
pw = data;
%%

for i = 1:6
P.v = pw.path_waves.aorta_finger(i).P{32};
P.fs = pw.path_waves.fs;
P_out.v = pw.config.p_out(i,:)*ones(size(P.v));
P_out.fs = pw.path_waves.fs  ;
Q.v = pw.path_waves.aorta_finger(i).U{32}.*pw.path_waves.aorta_finger(i).A{32};
Q.fs = pw.path_waves.fs;
curr_PPG_WK = estimate_ppg_using_windkessel(P, Q, P_out);
curr_PPG_WK.v = curr_PPG_WK.v(:);
[~,rel_el] = min(curr_PPG_WK.v);
curr_PPG_WK.v = [curr_PPG_WK.v(rel_el:end); 
curr_PPG_WK.v(1:rel_el-1)];
% find normalisation scaling factor
norm_factor = 1./range(curr_PPG_WK.v);
% normalise
curr_PPG_WK.v = (curr_PPG_WK.v-min(curr_PPG_WK.v)).*norm_factor;

data.waves.P_Radial{i} = pw.path_waves.aorta_finger(i).P{32};
data.waves.A_Radial{i} = pw.path_waves.aorta_finger(i).A{32};
data.waves.U_Radial{i} = pw.path_waves.aorta_finger(i).U{32};
data.waves.PPG_Radial{i} = curr_PPG_WK.v;
end
%%
% Setup figure for Pressure and Flow vel
paper_size = [350, 600];
figure('Position', [20,20,paper_size])
fig_settings.lwidth = 1.5;
fig_settings.ftsize = 14;
fig_settings.req_sites =  {'Carotid', 'AorticRoot', 'Brachial', 'Radial', 'Digital', 'Femoral', 'AntTibial'};
fig_settings.sig_types = {'P', 'U'};
fig_settings.colors = ["#0021F5","#8B8FF7","#CCCDFB","#F7CECD","#F0928F","#EB3323"];
fig_settings.legend = {'Age 25','Age 35','Age 45','Age 55','Age 65','Age 75'};
fig_settings.ylims = {[60,135], [-0.15, 0.75]};

% Plot
plot_baseline_for_age_signals(fig_settings, data)

% Save figure
PrintFigs(gcf, paper_size/70, PATHS.baseline_for_age_waves_fig_a)

% Setup figure for PPG and Area
figure('Position', [20,20,paper_size])
fig_settings.req_sites =  {'Carotid','AorticRoot','Brachial', 'Radial', 'Digital', 'Femoral', 'AntTibial'};
fig_settings.sig_types = {'A', 'PPG'};
fig_settings.colors = ["#0021F5","#8B8FF7","#CCCDFB","#F7CECD","#F0928F","#EB3323"];
fig_settings.legend = {'Age 25','Age 35','Age 45','Age 55','Age 65','Age 75'};
fig_settings.ylims = {'auto', [0, 1.1]};

% Plot
plot_baseline_for_age_signals(fig_settings, data)

% Save figure
PrintFigs(gcf, paper_size/70, PATHS.baseline_for_age_waves_fig_b)


%%
function plot_baseline_for_age_signals(fig_settings, data)

% identify baseline simulation data
baseline_sim_no = find(data.config.baseline_sim_for_age);

% cycle through different signals
for sig_type_no = 1 : length(fig_settings.sig_types)
    curr_sig_type = fig_settings.sig_types{sig_type_no};
    
    % cycle through different sites
    for req_site_no = 1 : length(fig_settings.req_sites)
        curr_site = fig_settings.req_sites{req_site_no};
        rel_sig.total_v = [];
        % extract relevant signal at this site
        for req_age_no = 1:length(baseline_sim_no)
            curr_age = baseline_sim_no(req_age_no);
            eval(['rel_sig.v = data.waves.' curr_sig_type '_' curr_site '{curr_age};'])
            rel_sig.fs = data.waves.fs;
            rel_sig.t = [0:length(rel_sig.v)-1]/rel_sig.fs;

            % convert to friendly units if needed
            if sum(strcmp(curr_sig_type, 'P'))
                rel_sig.units = 'mmHg';
                rel_sig.labels = 'Pressure';
                rel_sig.ylim = [60 90 120];
            elseif strcmp(curr_sig_type, 'A')
                rel_sig.v = sqrt(rel_sig.v*1000*1000/pi);   % m^2 to mm^2
                rel_sig.units = 'mm';
                rel_sig.labels = 'Luminal Diameter';
           
            elseif strcmp(curr_sig_type, 'PPG')
                rel_sig.units = 'au';
                rel_sig.labels = 'PPG';
                rel_sig.ylim = [0 0.5 1];
            elseif strcmp(curr_sig_type, 'U')
                rel_sig.units = 'm/s';
                rel_sig.labels = 'Flow Velocity';
                rel_sig.ylim = [0 0.3 0.6];
            end

            % setup subplot
            if sig_type_no == 1
                subplot(length(fig_settings.req_sites),2,(2*req_site_no)-1)
            else
                subplot(length(fig_settings.req_sites),2,(2*req_site_no))
            end
            rel_sig.total_v = [rel_sig.total_v ;rel_sig.v];
            % plot signal
            plot(rel_sig.t, rel_sig.v, 'Color', fig_settings.colors(curr_age), 'LineWidth', fig_settings.lwidth)
            hold on
        end
            % tidy up
            set(gca, 'XTick', 0:0.2:1)
            set(gca, 'FontSize', fig_settings.ftsize,'XAxisLocation', 'origin')
            xlim([0, rel_sig.t(end)])
            if req_site_no == 1
%                 dim = [.2+.5*(sig_type_no-1) .77 .2 .2];
                str = [rel_sig.labels, ' [', rel_sig.units, ']'];
%                 annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle', 'none', 'FontSize', fig_settings.ftsize);
                title(str, 'FontSize', fig_settings.ftsize);
            end
            if req_site_no < length(fig_settings.req_sites)
                set(gca, 'XTickLabel', [])
            end
            if sum(strcmp(fieldnames(fig_settings), 'ylims'))
                ylims = fig_settings.ylims{sig_type_no};
                if ~isnumeric(ylims)
                    temp_range = range(rel_sig.total_v);
                    ylims = [min(rel_sig.total_v)-0.2*temp_range, max(rel_sig.total_v)+0.2*temp_range];
                    ylim(ylims)
%                     step = round(max(rel_sig.total_v)-min(rel_sig.total_v)+0.4*temp_range,1);
                    if max(rel_sig.total_v) - min(rel_sig.total_v)+0.4*temp_range  > 3 && max(rel_sig.total_v) - min(rel_sig.total_v)+0.4*temp_range  < 100
                        yt1 = floor(min(rel_sig.total_v));
                        yt3 = ceil(max(rel_sig.total_v));
                        if mod(yt1,2) == 1
                            yt1 = yt1-1;
                        end
                        if mod(yt3,2) == 1
                            yt3 = yt3+1;
                        end
                        yticks([yt1, (yt1+yt3)/2 ,yt3])
                    elseif max(rel_sig.total_v) - min(rel_sig.total_v)+0.4*temp_range  <= 3
                        yt1 = round(min(rel_sig.total_v),1);
                        yt3 = round(max(rel_sig.total_v),1);

                        if mod(yt1*10,2) == 1
                            yt1 = (yt1*10-1)/10;
                        end
                        if mod(yt3*10,2) == 1
                            yt3 = (yt3*10-1)/10;
                        end
                        if yt1 < min(rel_sig.total_v)-0.2*temp_range
                            ylim([yt1-0.1*temp_range,max(rel_sig.total_v)+0.2*temp_range])
                        end
                        if yt3 > max(rel_sig.total_v)+0.2*temp_range
                            ylim([min(rel_sig.total_v)-0.2*temp_range,yt3+0.1*temp_range])
                        end
                        if yt3 == yt1
                           k =  (yt1+yt3)/2;
                           yt1 = yt1 - k;
                           yt3 = yt3 + k;                           
                        end
                        yticks([yt1, (yt1+yt3)/2 ,yt3])
                    else
                        yt1 = round(min(rel_sig.total_v)-0.1*temp_range,-1);
                        yt3 = round(max(rel_sig.total_v)+0.1*temp_range,-1);
                        yticks([yt1, (yt1+yt3)/2 ,yt3])
                    end
                else
                    ylim(ylims)
                    yticks(rel_sig.ylim)
                end
         end
            
            if req_site_no == length(fig_settings.req_sites)
                xlab = xlabel('Time [s]', 'FontSize', fig_settings.ftsize);
                xlim([0,0.9])
                xticks([0 0.2 0.4 0.6 0.8])
                xtickangle(45)
                if ylims(1)<0                    
                    set(xlab, 'Units', 'Normalized', 'Position', [0.75, -0.92, 0]); 
                end
            end
            box off
    end
end
end

%%
function PrintFigs(h, paper_size, savepath, close_fig)

if nargin<4
    close_fig = true;
end

set(h,'PaperUnits','inches');
set(h,'PaperSize', [paper_size(1), paper_size(2)]);
set(h,'PaperPosition',[0 0 paper_size(1) paper_size(2)]);
set(gcf,'color','w');
% print(h,'-dpdf',savepath)
% print(h,'-depsc',savepath)
print(h,'-dpng',savepath)

% if you want .eps illustrations, then do as follows:
up.eps_figs = 0;
if up.eps_figs
    % you need to download 'export_fig' from:
    % http://uk.mathworks.com/matlabcentral/fileexchange/23629-export-fig
    export_fig_dir_path = 'C:\Documents\Google Drive\Work\Projects\PhD\Github\phd\Tools\Other Scripts\export_fig\altmany-export_fig-76bd7fa\';
    addpath(export_fig_dir_path)
    export_fig(savepath, '-eps')
end
if close_fig
    close all;
end

% save 
fid = fopen([savepath, '.txt'], 'w');
p = mfilename('fullpath');
p = strrep(p, '\', '\\');
fprintf(fid, ['Figures generated by:\n\n ' p '.m \n\n on ' datestr(today)]);
fclose all;

end

%%
function ppg = estimate_ppg_using_windkessel(P, Q, Pout)

%% Calculate PPG

% Calculate time vectors for input signals
P.t = [0:length(P.v)-1]/P.fs;
Q.t = [0:length(Q.v)-1]/Q.fs;

% Find the resistance to flow further down the arterial tree
temp = P.v-Pout.v; % pressure drop between this segment and end of arterial tree
R = sum(temp)/sum(Q.v); % resistance is mean pressure drop over mean flow

% Find the flow into the more distal part of the arterial tree from this segment
Qout.v = (P.v-Pout.v)./R;  % I = V/R (electrical circuit)
Qout.t = P.t;

% Find the volume stored in the arterial segment
Volvb.t = Q.t;
const = 0;
Volvb.v = const + cumsum(Q.v) - cumsum(Qout.v);  % volume stored is the difference between inflow and outflow

ppg.t = Volvb.t;
% ppg.v = normwave(Volvb.v);
ppg.v = Volvb.v;

end
