% clear all;clc
% PATH_ROOT = 'C:\Users\43897\';
 PATH_ROOT = 'C:\Users\Jingyuan Hong\';
% PATH_ROOT = 'C:\Users\jh22\';
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
load([PATH_ROOT,'Downloads\pwdb_data_w_aorta_finger_path.mat'])
% PATH_ROOT = '/Users/jingyuanhong/';
% load([PATH_ROOT,'Downloads/pwdb_data_w_aorta_finger_path.mat']) % for mac

%% aoPWV
% extract aoPWV from flow velocity between ascending and descending aorta
curr_pwv = zeros(length(data.haemods),1);
for sim_no = 1: length(data.haemods)
    wave1.v = data.path_waves.aorta_foot(sim_no).U{1, 1};   % seg. 01-ascending aorta
    wave1.t = data.path_waves.aorta_foot(sim_no).onset_time(1);
    wave2.v = data.path_waves.aorta_foot(sim_no).U{1, 19};  % seg. 18-descending aorta
    wave2.t = data.path_waves.aorta_foot(sim_no).onset_time(19);
    lens = [length(wave1.v), length(wave2.v)];
    %%
        if abs(diff(lens))/min(lens) > 0.1
            fprintf(['Different length PWs for sim ' num2str(sim_no) ', ' curr_pwv_type])
            Foot_TT = nan;
        else
            % time-align waves
            % check waves are same duration
            if length(wave2.v) == length(wave1.v)+1
                wave2.v = wave2.v(1:end-1);
            elseif length(wave2.v) == length(wave1.v)+2
                wave2.v = wave2.v(1:end-2);
            elseif length(wave1.v) == length(wave2.v)+1
                wave1.v = wave1.v(1:end-1);
            elseif length(wave1.v) == length(wave2.v)+2
                wave1.v = wave1.v(1:end-2);
            end
            len = length(wave1.v);
            if length(wave2.v) ~= length(wave1.v)
                fprintf(['Different length PWs for sim ' num2str(sim_no) ', ' curr_pwv_type])
                Foot_TT = nan;
            else
                % repeat waves
                wave1.v = (wave1.v - min(wave1.v))/(max(wave1.v)-min(wave1.v)); %normalization
                wave2.v = (wave2.v - min(wave2.v))/(max(wave2.v)-min(wave2.v)); %normalization
                temp = linspace(wave1.v(1),wave1.v(end),length(wave1.v));
                wave1.v = wave1.v+ wave1.v(1)-temp(:);
                wave1.v = repmat(wave1.v, [5,1]);
                temp = linspace(wave2.v(1),wave2.v(end),length(wave2.v));
                wave2.v = wave2.v+wave2.v(1)-temp(:);
                wave2.v = repmat(wave2.v, [5,1]); 
            
            end
        end
Foot_TT = TTAlgorithm([wave1.v wave2.v],500,1,2,1,0);
lens_ascen_descen_aorta = data.config.network.length(sim_no,[1,2,14,18,27]);
path_len = sum(lens_ascen_descen_aorta);
curr_pwv(sim_no) = path_len./(Foot_TT(1)+wave2.t-wave1.t);
s = curr_pwv(sim_no);
pro = sprintf('Processing %i/%i\n',[sim_no,length(data.haemods)]);
fprintf(pro)
clc
end
aoPWV = curr_pwv;
%%
save aoPWV.mat aoPWV
clearvars -except data PATH_ROOT
%% cfPWV
curr_pwv = zeros(length(data.haemods),1);
for sim_no = 1 : length(data.haemods)
    wave1.v = data.path_waves.aorta_brain(sim_no).P{1, 10}; % seg. 15-carotid aorta-half
    wave1.t = data.path_waves.aorta_brain(sim_no).onset_time(10);
    wave2.v = data.path_waves.aorta_foot(sim_no).P{1, 55};  % seg. 46-femoral aorta-half
    wave2.t = data.path_waves.aorta_foot(sim_no).onset_time(55);
    lens = [length(wave1.v), length(wave2.v)];
    %%
        if abs(diff(lens))/min(lens) > 0.1
            fprintf(['Different length PWs for sim ' num2str(sim_no) ', ' curr_pwv_type])
            Foot_TT = nan;
        else     
            % time-align waves
            % check waves are same duration
            if length(wave2.v) == length(wave1.v)+1
                wave2.v = wave2.v(1:end-1);
            elseif length(wave2.v) == length(wave1.v)+2
                wave2.v = wave2.v(1:end-2);
            elseif length(wave1.v) == length(wave2.v)+1
                wave1.v = wave1.v(1:end-1);
            elseif length(wave1.v) == length(wave2.v)+2
                wave1.v = wave1.v(1:end-2);
            end
            if length(wave2.v) ~= length(wave1.v)
                fprintf(['Different length PWs for sim ' num2str(sim_no) ', ' curr_pwv_type])
                Foot_TT = nan;
            else
                % repeat waves
                wave1.v = (wave1.v - min(wave1.v))/(max(wave1.v)-min(wave1.v)); %normalization
                wave2.v = (wave2.v - min(wave2.v))/(max(wave2.v)-min(wave2.v)); %normalization
                temp = linspace(wave1.v(1),wave1.v(end),length(wave1.v));
                wave1.v = wave1.v+ wave1.v(1)-temp(:);
                wave1.v = repmat(wave1.v, [5,1]);
                temp = linspace(wave2.v(1),wave2.v(end),length(wave2.v));
                wave2.v = wave2.v+wave2.v(1)-temp(:);
                wave2.v = repmat(wave2.v, [5,1]); 

            end
        end
Foot_TT = TTAlgorithm([wave1.v wave2.v],500,1,1,1,0); % sim_no = 1677 should use algorithm 1 to avoid large error % close plot function when you are compiling

lens_aorta_femoral = [data.config.network.length(sim_no,[1,2,14,18,27,28,35,37,39,41,42,44]), data.config.network.length(sim_no,46)/2];
lens_aorta_carotid = [data.config.network.length(sim_no,[1,2]), data.config.network.length(sim_no,15)/2];
path_len = sum(lens_aorta_femoral)-sum(lens_aorta_carotid);

curr_pwv(sim_no) = path_len./(Foot_TT(1)+wave2.t-wave1.t);
end
cfPWV = curr_pwv;
save cfPWV.mat cfPWV
clearvars -except data PATH_ROOT
%% baPWV
curr_pwv = zeros(length(data.haemods),1);
for sim_no = 1 : length(data.haemods)
    wave1.v = data.path_waves.aorta_finger(sim_no).P{1, 24};  % seg. 21-brachial artery
    wave1.t = data.path_waves.aorta_finger(sim_no).onset_time(24);
    wave2.v = data.path_waves.aorta_foot(sim_no).P{1, 72};  % seg. 49-ankle/anterior tibial artery
    wave2.t = data.path_waves.aorta_foot(sim_no).onset_time(72);
    lens = [length(wave1.v), length(wave2.v)];
    %%
        if abs(diff(lens))/min(lens) > 0.1
            fprintf(['Different length PWs for sim ' num2str(sim_no) ', ' curr_pwv_type])
            Foot_TT_ba = nan;
        else
            
            % time-align waves
            % check waves are same duration
            if length(wave2.v) == length(wave1.v)+1
                wave2.v = wave2.v(1:end-1);
            elseif length(wave2.v) == length(wave1.v)+2
                wave2.v = wave2.v(1:end-2);
            elseif length(wave1.v) == length(wave2.v)+1
                wave1.v = wave1.v(1:end-1);
            elseif length(wave1.v) == length(wave2.v)+2
                wave1.v = wave1.v(1:end-2);
            end
            if length(wave2.v) ~= length(wave1.v)
                fprintf(['Different length PWs for sim ' num2str(sim_no) ', ' curr_pwv_type])
                Foot_TT_ba = nan;
            else
%                 % repeat waves
                wave1.v = (wave1.v - min(wave1.v))/(max(wave1.v)-min(wave1.v)); %normalization
                wave2.v = (wave2.v - min(wave2.v))/(max(wave2.v)-min(wave2.v)); %normalization
                temp = linspace(wave1.v(1),wave1.v(end),length(wave1.v));
                wave1.v = wave1.v+ wave1.v(1)-temp(:);
                wave1.v = repmat(wave1.v, [5,1]);
                temp = linspace(wave2.v(1),wave2.v(end),length(wave2.v));
                wave2.v = wave2.v+wave2.v(1)-temp(:);
                wave2.v = repmat(wave2.v, [5,1]); 

            end
        end
Foot_TT_ba = TTAlgorithm([wave1.v wave2.v],500,1,1,1,0);
    
lens_aortic_brachial = [data.config.network.length(sim_no,[1,2,14,19]),data.config.network.length(sim_no,21)*0.75];
lens_aortic_ankle = data.config.network.length(sim_no,[1,2,14,18,27,28,35,37,39,41,42,44,46,49]);

path_len = sum(lens_aortic_ankle) - sum(lens_aortic_brachial);
curr_pwv(sim_no) = path_len./(Foot_TT_ba(1)+wave2.t-wave1.t);
end
baPWV = curr_pwv;
save baPWV.mat baPWV
clearvars -except data PATH_ROOT
%% cbPWV
curr_pwv = zeros(length(data.haemods),1);
for sim_no = 1 : length(data.haemods)
    wave1.v = data.path_waves.aorta_brain(sim_no).P{1, 10};   % seg. 15-carotid aorta-half
    wave1.t = data.path_waves.aorta_brain(sim_no).onset_time(10);
    wave2.v = data.path_waves.aorta_finger(sim_no).P{1, 24};  % seg. 21-brachial aorta-3/4
    wave2.t = data.path_waves.aorta_finger(sim_no).onset_time(24);
    lens = [length(wave1.v), length(wave2.v)];
    %%
        if abs(diff(lens))/min(lens) > 0.1
            fprintf(['Different length PWs for sim ' num2str(sim_no) ', ' curr_pwv_type])
            Foot_TT = nan;
        else
            
            % time-align waves
            % check waves are same duration
            if length(wave2.v) == length(wave1.v)+1
                wave2.v = wave2.v(1:end-1);
            elseif length(wave2.v) == length(wave1.v)+2
                wave2.v = wave2.v(1:end-2);
            elseif length(wave1.v) == length(wave2.v)+1
                wave1.v = wave1.v(1:end-1);
            elseif length(wave1.v) == length(wave2.v)+2
                wave1.v = wave1.v(1:end-2);
            end
            if length(wave2.v) ~= length(wave1.v)
                fprintf(['Different length PWs for sim ' num2str(sim_no) ', ' curr_pwv_type])
                Foot_TT = nan;
            else
                % repeat waves
                wave1.v = (wave1.v - min(wave1.v))/(max(wave1.v)-min(wave1.v)); %normalization
                wave2.v = (wave2.v - min(wave2.v))/(max(wave2.v)-min(wave2.v)); %normalization
                temp = linspace(wave1.v(1),wave1.v(end),length(wave1.v));
                wave1.v = wave1.v+ wave1.v(1)-temp(:);
                wave1.v = repmat(wave1.v, [5,1]);
                temp = linspace(wave2.v(1),wave2.v(end),length(wave2.v));
                wave2.v = wave2.v+wave2.v(1)-temp(:);
                wave2.v = repmat(wave2.v, [5,1]); 

            end
        end
Foot_TT = TTAlgorithm([wave1.v wave2.v],500,1,1,1,0);
lens_aortic_brachial = [data.config.network.length(sim_no,[1,2,14,19]),data.config.network.length(sim_no,21)*0.75];
lens_aorta_carotid = [data.config.network.length(sim_no,[1,2]), data.config.network.length(sim_no,15)/2];
path_len = sum(lens_aortic_brachial)-sum(lens_aorta_carotid);

curr_pwv(sim_no) = path_len./(Foot_TT(1)+wave2.t-wave1.t);
end
cbPWV = curr_pwv;
save cbPWV.mat cbPWV
clearvars -except data PATH_ROOT
%% crPWV
curr_pwv = zeros(length(data.haemods),1);
for sim_no = 1 : length(data.haemods)
    wave1.v = data.path_waves.aorta_brain(sim_no).P{1, 10};   % seg. 15-carotid aorta-half
    wave1.t = data.path_waves.aorta_brain(sim_no).onset_time(10);
    wave2.v = data.path_waves.aorta_finger(sim_no).P{1, 32};  % seg. 22-radial aorta
    wave2.t = data.path_waves.aorta_finger(sim_no).onset_time(32);
    lens = [length(wave1.v), length(wave2.v)];
    %%
        if abs(diff(lens))/min(lens) > 0.1
            fprintf(['Different length PWs for sim ' num2str(sim_no) ', ' curr_pwv_type])
            Foot_TT = nan;
        else           
            % time-align waves
            % check waves are same duration
            if length(wave2.v) == length(wave1.v)+1
                wave2.v = wave2.v(1:end-1);
            elseif length(wave2.v) == length(wave1.v)+2
                wave2.v = wave2.v(1:end-2);
            elseif length(wave1.v) == length(wave2.v)+1
                wave1.v = wave1.v(1:end-1);
            elseif length(wave1.v) == length(wave2.v)+2
                wave1.v = wave1.v(1:end-2);
            end
            if length(wave2.v) ~= length(wave1.v)
                fprintf(['Different length PWs for sim ' num2str(sim_no) ', ' curr_pwv_type])
                Foot_TT = nan;
            else
                % repeat waves
                wave1.v = (wave1.v - min(wave1.v))/(max(wave1.v)-min(wave1.v)); %normalization
                wave2.v = (wave2.v - min(wave2.v))/(max(wave2.v)-min(wave2.v)); %normalization
                temp = linspace(wave1.v(1),wave1.v(end),length(wave1.v));
                wave1.v = wave1.v+ wave1.v(1)-temp(:);
                wave1.v = repmat(wave1.v, [5,1]);
                temp = linspace(wave2.v(1),wave2.v(end),length(wave2.v));
                wave2.v = wave2.v+wave2.v(1)-temp(:);
                wave2.v = repmat(wave2.v, [5,1]); 

            end
        end
Foot_TT = TTAlgorithm([wave1.v wave2.v],500,1,1,1,0);
lens_aorta_radial = [data.config.network.length(sim_no,[1,2,14,19,21]),data.config.network.length(sim_no,22)/2];
lens_aorta_carotid = [data.config.network.length(sim_no,[1,2]), data.config.network.length(sim_no,15)/2];
path_len = sum(lens_aorta_radial)-sum(lens_aorta_carotid);

curr_pwv(sim_no) = path_len./(Foot_TT(1)+wave2.t-wave1.t);
end
crPWV = curr_pwv;
save crPWV.mat crPWV
clearvars -except data PATH_ROOT
%% ftPWV
curr_pwv = zeros(length(data.haemods),1);
for sim_no = 1 : length(data.haemods)
    
    % P.v = data.path_waves.aorta_finger(sim_no).P{1,41};% seg. 112-index finger artery
    % P.fs = data.path_waves.fs;
    % P_out.v = data.config.p_out(sim_no,:)*ones(size(P.v));
    % P_out.fs = data.path_waves.fs  ;
    % Q.v = data.path_waves.aorta_finger(sim_no).U{41}.*data.path_waves.aorta_finger(sim_no).A{41};
    % Q.fs = data.path_waves.fs;
    % curr_PPG_WK = estimate_ppg_using_windkessel(P, Q, P_out);
    % curr_PPG_WK.v = curr_PPG_WK.v(:);
    % [~,rel_el] = min(curr_PPG_WK.v);
    % curr_PPG_WK.v = [curr_PPG_WK.v(rel_el:end); 
    % curr_PPG_WK.v(1:rel_el-1)];
    % % find normalisation scaling factor
    % norm_factor = 1./range(curr_PPG_WK.v);
    % % normalise
    % wave1.v = (curr_PPG_WK.v-min(curr_PPG_WK.v)).*norm_factor;
    % clear P P_out Q curr_PPG_WK norm_factor
    % 
    % wave1.t = data.path_waves.aorta_finger(sim_no).onset_time(41);
    load([PATH_ROOT,'OneDrive - King''s College London\PhD Work\code\PW_studies\pwdb_1\PWs_mat\mat\PWs_Digital.mat']);
    wave1.v = PWs.PPG{sim_no};
    wave1.t = PWs.onset_times.PPG(sim_no);

    % P.v = data.path_waves.aorta_foot(sim_no).P{1,72}; %seg. 49-ankle/anterior tibial artery
    % P.fs = data.path_waves.fs;
    % P_out.v = data.config.p_out(sim_no,:)*ones(size(P.v));
    % P_out.fs = data.path_waves.fs  ;
    % Q.v = data.path_waves.aorta_foot(sim_no).U{72}.*data.path_waves.aorta_foot(sim_no).A{72};
    % Q.fs = data.path_waves.fs;
    % curr_PPG_WK = estimate_ppg_using_windkessel(P, Q, P_out);
    % curr_PPG_WK.v = curr_PPG_WK.v(:);
    % [~,rel_el] = min(curr_PPG_WK.v);
    % curr_PPG_WK.v = [curr_PPG_WK.v(rel_el:end); 
    % curr_PPG_WK.v(1:rel_el-1)];
    % % find normalisation scaling factor
    % norm_factor = 1./range(curr_PPG_WK.v);
    % % normalise
    % wave2.v = (curr_PPG_WK.v-min(curr_PPG_WK.v)).*norm_factor;
    % 
    % clear P P_out Q curr_PPG_WK norm_factor
    % wave2.t = data.path_waves.aorta_foot(sim_no).onset_time(72);

    load([PATH_ROOT,'OneDrive - King''s College London\PhD Work\code\PW_studies\pwdb_1\PWs_mat\mat\PWs_AntTibial.mat']);
    wave2.v = PWs.PPG{sim_no};
    wave1.t = PWs.onset_times.PPG(sim_no);
    lens = [length(wave1.v), length(wave2.v)];
    %%
        if abs(diff(lens))/min(lens) > 0.1
            fprintf(['Different length PWs for sim ' num2str(sim_no)])
            Foot_TT_ba = nan;
        else
            
            % time-align waves
            % check waves are same duration
            if length(wave2.v) == length(wave1.v)+1
                wave2.v = wave2.v(1:end-1);
            elseif length(wave2.v) == length(wave1.v)+2
                wave2.v = wave2.v(1:end-2);
            elseif length(wave1.v) == length(wave2.v)+1
                wave1.v = wave1.v(1:end-1);
            elseif length(wave1.v) == length(wave2.v)+2
                wave1.v = wave1.v(1:end-2);
            end
            if length(wave2.v) ~= length(wave1.v)
                fprintf(['Different length PWs for sim ' num2str(sim_no) ', ' curr_pwv_type])
                Foot_TT_ba = nan;
            else
%                 % repeat waves
                wave1.v = (wave1.v - min(wave1.v))/(max(wave1.v)-min(wave1.v)); %normalization
                wave2.v = (wave2.v - min(wave2.v))/(max(wave2.v)-min(wave2.v)); %normalization
                temp = linspace(wave1.v(1),wave1.v(end),length(wave1.v));
                wave1.v = wave1.v+ wave1.v(1)-temp(:);
                wave1.v = repmat(wave1.v, [5,1]);
                temp = linspace(wave2.v(1),wave2.v(end),length(wave2.v));
                wave2.v = wave2.v+wave2.v(1)-temp(:);
                wave2.v = repmat(wave2.v, [5,1]); 

            end
        end
Foot_TT_da = TTAlgorithm([wave1.v wave2.v],500,1,1,1,0);% 2 is better to extract TT
    
lens_aortic_digital = data.config.network.length(sim_no,[1,2,14,19,21,22,108,112]);
lens_aortic_ankle = data.config.network.length(sim_no,[1,2,14,18,27,28,35,37,39,41,42,44,46,49]);

path_len = sum(lens_aortic_ankle) - sum(lens_aortic_digital);

curr_pwv(sim_no) = path_len./(Foot_TT_da(1)+wave2.t-wave1.t);
end
ftPWV = curr_pwv;
save ftPWV.mat ftPWV
clearvars -except data PATH_ROOT
%% cPP bPP
cPP = zeros(length(data.haemods),1);
bPP = zeros(length(data.haemods),1);
for sim_no = 1:length(data.haemods)
P_bra = data.path_waves.aorta_finger(sim_no).P{1, 24};
P_aor = data.path_waves.aorta_finger(sim_no).P{1, 1};
cPP(sim_no) = max(P_aor) - min(P_aor);
bPP(sim_no) = max(P_bra) - min(P_bra);
end
save PP.mat cPP bPP
clearvars -except data PATH_ROOT
%% AP AIx
options.do_plot = 0;
options.exclude_low_quality_data = 1;
AIx = zeros(length(data.haemods),1);
AP = zeros(length(data.haemods),1);
for sim_no = 1:length(data.haemods)
    S.v = data.path_waves.aorta_finger(sim_no).P{1, 1};
    S.fs = 500;
    [~, fid_pts, ~, ~] = PulseAnalyse10(S, options);
    p1in = fid_pts.p1in;
    p2pk = fid_pts.p2pk;
    if isnan(p1in)||isnan(p2pk)
        fprintf('Wrong!\n')
    else
        AP(sim_no) = S.v(p2pk) - S.v(p1in);
        PP = max(S.v) - min(S.v);
        AIx(sim_no) = AP(sim_no)/PP*100;
    end
end
save AIx_AP.mat AP AIx
clearvars -except data PATH_ROOT
%% Pb RM
wave_a.Pf_A = zeros(length(data.haemods),1);
wave_a.Pb_A = wave_a.Pf_A;
wave_a.RM = wave_a.Pf_A;
for sim_no = 1 : length(data.haemods)

    wave_a.p = data.path_waves.aorta_finger(sim_no).P{1, 1};  % seg. 1-aorta
    wave_a.v = data.path_waves.aorta_finger(sim_no).U{1, 1}; % seg. 1-aorta

    wave_a.rho = data.config.constants.rho(sim_no);
    wave_a.c = sqrt(sum(diff(133.322*wave_a.p).^2)./sum(diff(wave_a.v).^2))/wave_a.rho;

    wave_a.dP_f = (diff(133.322*wave_a.p) + (wave_a.rho*wave_a.c.*diff(wave_a.v)))/2;
    wave_a.dP_b = (diff(133.322*wave_a.p) - (wave_a.rho.*wave_a.c.*diff(wave_a.v)))/2;
    wave_a.P_f = cumsum(wave_a.dP_f)/133.322;
    wave_a.P_b = cumsum(wave_a.dP_b)/133.322;
    wave_a.Pf_A(sim_no) = max(wave_a.P_f);
    wave_a.Pb_A(sim_no) = max(wave_a.P_b);
    wave_a.RM(sim_no) = wave_a.Pb_A(sim_no)/wave_a.Pf_A(sim_no);
end
Pb = wave_a.Pb_A;
RM = wave_a.RM;
save Pb_RM.mat Pb RM
clearvars -except data PATH_ROOT
%% RIppg
load([PATH_ROOT,'OneDrive - King''s College London\PhD Work\code\PW_studies\pwdb_1\PWs_mat\mat\PWs_Digital.mat'])
PWs_indfinger = PWs;
diaT_PPG = zeros(length(PWs_indfinger.PPG),1);
dia_PPG = zeros(length(PWs_indfinger.PPG),1);
max_PPG = dia_PPG;
RI_PPG = dia_PPG;
for i = 1: length(PWs_indfinger.PPG)
    sig.v = PWs_indfinger.PPG{1, i};
    sig.fs = 500;
    options.do_plot = 0; sig.ht = 1.75;
    [cv_inds, ~, ~, ~] = PulseAnalyse10(sig, options);
    RI_PPG(i) = cv_inds.RI.v;

end
save RI_PPG.mat RI_PPG
clearvars -except data PATH_ROOT
%% CAVI
curr_ind = zeros(length(data.haemods),1);
ha_pwv = curr_ind; % heart-to-ankle PWV
ha_beta = curr_ind;
a = curr_ind;
b = curr_ind;
ind_1 = [0.850 0.658 0.432];
ind_2 = [0.695 2.103 4.441];

Ps = zeros(length(data.haemods),1);
Pd = zeros(length(data.haemods),1);
Pp = zeros(length(data.haemods),1);
for sim_no = 1:length(data.haemods)
P_bra = data.path_waves.aorta_finger(sim_no).P{1, 24};
Ps(sim_no) = max(P_bra);
Pd(sim_no) = min(P_bra);
Pp(sim_no) = max(P_bra) - min(P_bra);
end

for sim_no = 1 : length(data.haemods)
    wave1.v = data.path_waves.aorta_foot(sim_no).P{1, 1};   % seg. 01-aorta
    wave1.t = data.path_waves.aorta_foot(sim_no).onset_time(1);
    wave2.v = data.path_waves.aorta_foot(sim_no).P{1, 72};  % seg. 49-ankle/anterior tibial artery
    wave2.t = data.path_waves.aorta_foot(sim_no).onset_time(72);
    lens_ab = [length(wave1.v), length(wave2.v)];
    
    %%
        if abs(diff(lens_ab))/min(lens_ab) > 0.1
            fprintf(['Different length PWs for sim ' num2str(sim_no) ', ' curr_pwv_type])
            Foot_TT_ha = nan;
        else
            
            % time-align waves
            % check waves are same duration
            if length(wave2.v) == length(wave1.v)+1
                wave2.v = wave2.v(1:end-1);
            elseif length(wave2.v) == length(wave1.v)+2
                wave2.v = wave2.v(1:end-2);
            elseif length(wave1.v) == length(wave2.v)+1
                wave1.v = wave1.v(1:end-1);
            elseif length(wave1.v) == length(wave2.v)+2
                wave1.v = wave1.v(1:end-2);
            end
            if length(wave2.v) ~= length(wave1.v)
                fprintf(['Different length PWs for sim ' num2str(sim_no) ', ' curr_pwv_type])
                Foot_TT_ha = nan;
            else
%                 % repeat waves
                wave1.v = (wave1.v - min(wave1.v))/(max(wave1.v)-min(wave1.v)); %normalization
                wave2.v = (wave2.v - min(wave2.v))/(max(wave2.v)-min(wave2.v)); %normalization
                temp = linspace(wave1.v(1),wave1.v(end),length(wave1.v));
                wave1.v = wave1.v+ wave1.v(1)-temp(:);
                wave1.v = repmat(wave1.v, [5,1]);
                temp = linspace(wave2.v(1),wave2.v(end),length(wave2.v));
                wave2.v = wave2.v+wave2.v(1)-temp(:);
                wave2.v = repmat(wave2.v, [5,1]); 

            end
        end
Foot_TT_ha = TTAlgorithm([wave1.v wave2.v],500,1,1,1,0);
lens_aortic_ankle = data.config.network.length(sim_no,[1,2,14,18,27,28,35,37,39,41,42,44,46,49]);
path_len = sum(lens_aortic_ankle);
ha_pwv(sim_no) = path_len./(Foot_TT_ha(1)+wave2.t-wave1.t);
ha_beta(sim_no) = log(Ps(sim_no)/Pd(sim_no))*2*data.config.constants.rho(sim_no)*ha_pwv(sim_no)^2/(Pp(sim_no)*133.3224);% ln(mmHg*mmHg-1)*(kg*m-3)*(m2*s-2)/mmHg
% mmHg = 133.3224 Pa = 133.3224 kg*m-1*s-2

if ha_beta(sim_no) < 7.34875
    a(sim_no) = ind_1(1);
    b(sim_no) = ind_2(1);
elseif ha_beta(sim_no) >= 7.34875 && ha_beta(sim_no) < 10.30372
    a(sim_no) = ind_1(2);
    b(sim_no) = ind_2(2);
else
    a(sim_no) = ind_1(3);
    b(sim_no) = ind_2(3);
end

curr_ind(sim_no) = a(sim_no)*2*data.config.constants.rho(sim_no)*log(Ps(sim_no)/Pd(sim_no))*ha_pwv(sim_no)^2/(Pp(sim_no)*133.3224) + b(sim_no);
end
CAVI = curr_ind;
save CAVI.mat CAVI
clearvars -except data PATH_ROOT
%% DC
DC = zeros(length(data.haemods),1);
for sim_no = 1 : length(data.haemods)
    LA = data.path_waves.aorta_brain(sim_no).A{1, 10}; %carotid 10 
    PW = data.path_waves.aorta_finger(sim_no).P{1, 24};%brachial 24
    LAd = LA(1);
    LAs = max(LA);
    SBP = max(PW);
    DBP = PW(1);
    PP = (SBP - DBP)*0.13332;
    DC(sim_no) = 1e3*(LAs - LAd)/(LAd*PP);
end
save DC.mat DC
clearvars -except data PATH_ROOT
%% SI

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