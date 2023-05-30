clear all;clc
% PATH_ROOT = 'C:\Users\43897\';
PATH_ROOT = 'C:\Users\Jingyuan Hong\';
% PATH_ROOT = 'C:\Users\jh22\';
load([PATH_ROOT,'\Downloads\pwdb_data_w_aorta_finger_path'])
%%
load([PATH_ROOT,'OneDrive - King''s College London\PhD Work\code\PW_studies\pwdb_1\PWs_mat\mat\PWs_Digital.mat'])

%% Aortic Backward wave amplitude and reflection magnitude
wave_a.Pf_A = zeros(length(data.haemods),1);
wave_a.Pb_A = zeros(length(data.haemods),1);
wave_a.RM = zeros(length(data.haemods),1);
for sim_no = 1
    
    wave_a.p = data.path_waves.aorta_finger(sim_no).P{1, 1};  % seg. 1-aorta
    wave_a.v = data.path_waves.aorta_finger(sim_no).U{1, 1}; % seg. 1-aorta
   
    wave_a.rho = data.config.constants.rho(sim_no);
    wave_a.c = sqrt(sum(diff(133.322*wave_a.p).^2)./sum(diff(wave_a.v).^2))/wave_a.rho;
    
    wave_a.dP_f = (diff(133.322*wave_a.p) + (wave_a.rho*wave_a.c*diff(wave_a.v)))/2;
    wave_a.dP_b = (diff(133.322*wave_a.p) - (wave_a.rho*wave_a.c*diff(wave_a.v)))/2;
    wave_a.P_f = cumsum(wave_a.dP_f)/133.322;%+0.5*wave_a.p(1);
    wave_a.P_b = cumsum(wave_a.dP_b)/133.322;%+0.5*wave_a.p(1);
    wave_a.Pf_A(sim_no) = max(wave_a.P_f);
    wave_a.Pb_A(sim_no) = max(wave_a.P_b);
    wave_a.RM(sim_no) = wave_a.Pb_A(sim_no)/wave_a.Pf_A(sim_no);
    
end

%%
PWs_indfinger = PWs;
diaT_PPG = zeros(length(PWs_indfinger.PPG),1);
dia_PPG = zeros(length(PWs_indfinger.PPG),1);
max_PPG = dia_PPG;
RI_PPG = dia_PPG;
% RI_PPG_p = dia_PPG;
RI_PPG_p1pk = dia_PPG;
for i = 1%: length(PWs_indfinger.PPG)
    PPG = PWs_indfinger.PPG{1, i};
    dPPG = PPG(2:end) - PPG(1:end-1);
    sdPPG = dPPG(2:end) - dPPG(1:end-1);
    zeroPPG = [];
    for j = 1:length(sdPPG)-50
        if sdPPG(j) >= 0 && sdPPG(j+1)<0
            zeroPPG = [zeroPPG;j];
        end
    end
    diaT_PPG(i) = zeroPPG(end);
    dia_PPG(i) = PPG(diaT_PPG(i));
    max_PPG(i) = max(PPG);
    RI_PPG(i) = (max_PPG(i) - dia_PPG(i))/max_PPG(i);
    RI_PPG(i) = dia_PPG(i)/max_PPG(i);
%     RI_PPG_p1pk(i) = (p1pk(i) - dia_PPG(i))/p1pk(i);
end
%%
paper_size = [600, 300];
figure('Position', [20,20,paper_size])
fig_settings.lwidth = 1;
fig_settings.ftsize = 10;

subplot(1,2,1)
plot([1:length(wave_a.p)]/500,wave_a.p,'Color',0*[1 1 1],'LineWidth', fig_settings.lwidth)
hold on
plot([1:length(wave_a.P_f)]/500,wave_a.P_f,'--','Color',0*[1 1 1],'LineWidth', fig_settings.lwidth+0.2)
hold on
plot([1:length(wave_a.P_b)]/500,wave_a.P_b,'-.','Color',0*[1 1 1],'LineWidth', fig_settings.lwidth+0.2)
hold on
plot((find(wave_a.p == struct2table(data.haemods).P1in_a(1)))/500,struct2table(data.haemods).P1in_a(1)*1.02,'rv');%,'MarkerSize',8,'Linewidth',1)
hold on
plot((find(wave_a.p == struct2table(data.haemods).P2pk_a(1)))/500,struct2table(data.haemods).P2pk_a(1)*1.02,'rv');%,'MarkerSize',8,'Linewidth',1)
hold on
plot((find(wave_a.P_b == wave_a.Pb_A(1)))/500,wave_a.Pb_A(1)*1.15,'rv');%,'MarkerSize',8,'Linewidth',1)
hold on
plot((find(wave_a.P_f == wave_a.Pf_A(1)))/500,wave_a.Pf_A(1)*1.1,'rv');%,'MarkerSize',8,'Linewidth',1)
l1 = legend(gca,{'Central','Forward','Backward'},'FontSize',fig_settings.ftsize-2);
l1.Position(2) = 0.72;
l1.Position(1) = 0.30;


box off
xlabel('Time [s]','FontSize', fig_settings.ftsize-1)
ylabel('Pressure [mmHg]','FontSize', fig_settings.ftsize-1)

xlim([0,length(wave_a.p)/500])
ylim([-5,110])
set(gca, 'FontSize', fig_settings.ftsize-1)
breakyaxis([25 70]);

subplot(1,2,2)
PPG = PWs_indfinger.PPG{1, 1};
plot([1:length(PPG)]/500,PPG,'Color',0*[1 1 1],'LineWidth', fig_settings.lwidth)

xlim([0,length(sdPPG)/500])
ylim([0,1.05])
set(gca, 'FontSize', fig_settings.ftsize-1)
ylabel('PPG [au]','FontSize', fig_settings.ftsize-1)
xlabel('Time [s]','FontSize', fig_settings.ftsize-1)

box off
PATH_SAVE = [PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\figure\fig 3\'];
% saveas(gcf, [PATH_SAVE 'indices_plot_new.png'])