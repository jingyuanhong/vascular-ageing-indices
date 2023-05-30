for time = 1:5
for SNR = [15 20 30]
%% cfPWV Pressure 
fprintf('Add noise on cfPWV\n')
curr_pwv = zeros(length(data.haemods),1);
for sim_no = 1: length(data.haemods)
    rng(sim_no)
    seed = round(rand(1,2)*10000+time*SNR);
    wave1.v = noisegen(data.path_waves.aorta_brain(sim_no).P{1, 10},SNR,seed(1));   % seg. 15-carotid aorta-half
    wave1.t = data.path_waves.aorta_brain(sim_no).onset_time(10);
    wave2.v = noisegen(data.path_waves.aorta_foot(sim_no).P{1, 55},SNR,seed(2));    % seg. 46-femoral aorta-half
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
%%
filename = sprintf('cfPWV_n%i_%i.mat',[SNR,time]);
save(filename,'cfPWV')

%% aoPWV - Flow rate waveform
fprintf('Add noise on aoPWV\n')
curr_pwv = zeros(length(data.haemods),1);
for sim_no = 1 : length(data.haemods)
    rng(sim_no)
    seed = round(rand(1,2)*10000+time*SNR);
    wave1.v = noisegen(data.path_waves.aorta_foot(sim_no).U{1, 1},SNR,seed(1));   % seg. 01-ascending aorta
    wave1.t = data.path_waves.aorta_foot(sim_no).onset_time(1);
    wave2.v = noisegen(data.path_waves.aorta_foot(sim_no).U{1, 19},SNR,seed(2));  % seg. 27-descending aorta
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
lens_ascen_descen_aorta = data.config.network.length(sim_no,[1,2,14,18]);
path_len = sum(lens_ascen_descen_aorta);
curr_pwv(sim_no) = path_len./(Foot_TT(1)+wave2.t-wave1.t);
end

aoPWV = curr_pwv;
%%
filename = sprintf('aoPWV_n%i_%i.mat',[SNR,time]);
save(filename,'aoPWV')
%% Brachial-Ankle - Pressure or Volume
fprintf('Add noise on baPWV\n')
curr_pwv = zeros(length(data.haemods),1);
for sim_no = 1 : length(data.haemods)
    rng(sim_no)
    seed = round(rand(1,2)*10000+time*SNR);
    wave1.v = noisegen(data.path_waves.aorta_finger(sim_no).P{1, 24},SNR,seed(1));  % seg. 21-brachial artery
    wave1.t = data.path_waves.aorta_finger(sim_no).onset_time(24);
    wave2.v = noisegen(data.path_waves.aorta_foot(sim_no).P{1, 72},SNR,seed(2));  % seg. 49-ankle/anterior tibial artery
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
%%
filename = sprintf('baPWV_n%i_%i.mat',[SNR,time]);
save(filename,'baPWV')
%% Pressure crPWV
fprintf('Add noise on crPWV\n')
curr_pwv = zeros(length(data.haemods),1);
for sim_no = 1 : length(data.haemods)
    rng(sim_no)
    seed = round(rand(1,2)*10000+time*SNR);
    wave1.v = noisegen(data.path_waves.aorta_brain(sim_no).P{1, 10},SNR,seed(1));   % seg. 15-carotid aorta-half
    wave1.t = data.path_waves.aorta_brain(sim_no).onset_time(10);
    wave2.v = noisegen(data.path_waves.aorta_finger(sim_no).P{1, 32},SNR,seed(2));  % seg. 22-radial aorta
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
Foot_TT = TTAlgorithm([wave1.v wave2.v],500,1,1,1,0);% sim_no = 1677 should use algorithm 1 to avoid large error
% lens_aorta_radial = data.config.network.length(sim_no,[1,2,14,19,21,22]);
lens_aorta_radial = [data.config.network.length(sim_no,[1,2,14,19,21]),data.config.network.length(sim_no,22)/2];
lens_aorta_carotid = [data.config.network.length(sim_no,[1,2]), data.config.network.length(sim_no,15)/2];
path_len = sum(lens_aorta_radial)-sum(lens_aorta_carotid);

curr_pwv(sim_no) = path_len./(Foot_TT(1)+wave2.t-wave1.t);
end
crPWV = curr_pwv;

fprintf('Add noise on cbPWV\n')
curr_pwv = zeros(length(data.haemods),1);
for sim_no = 1 : length(data.haemods)
    rng(sim_no)
    seed = round(rand(1,2)*10000+time*SNR);
    wave1.v = noisegen(data.path_waves.aorta_brain(sim_no).P{1, 10},SNR,seed(1));   % seg. 15-carotid aorta-half
    wave1.t = data.path_waves.aorta_brain(sim_no).onset_time(10);
    wave2.v = noisegen(data.path_waves.aorta_finger(sim_no).P{1, 24},SNR,seed(2));  % seg. 21-brachial aorta-3/4
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
%%
filename = sprintf('cbcrPWV_n%i_%i.mat',[SNR,time]);
save(filename,'cbPWV','crPWV')
%% ftPWV
fprintf('Add noise on ftPWV\n')
curr_pwv = zeros(length(data.haemods),1);
for sim_no = 1 : length(data.haemods)
 % seg. 112-index finger artery
 % seg. 49-ankle/anterior tibial artery
    rng(sim_no)
    seed = round(rand(1,2)*10000+time*SNR);
    P.v = data.path_waves.aorta_finger(sim_no).P{1,41};
    P.fs = data.path_waves.fs;
    P_out.v = data.config.p_out(sim_no,:)*ones(size(P.v));
    P_out.fs = data.path_waves.fs  ;
    Q.v = data.path_waves.aorta_finger(sim_no).U{41}.*data.path_waves.aorta_finger(sim_no).A{41};
    Q.fs = data.path_waves.fs;
    curr_PPG_WK = estimate_ppg_using_windkessel(P, Q, P_out);
    curr_PPG_WK.v = curr_PPG_WK.v(:);
    [~,rel_el] = min(curr_PPG_WK.v);
    curr_PPG_WK.v = [curr_PPG_WK.v(rel_el:end); 
    curr_PPG_WK.v(1:rel_el-1)];
    % find normalisation scaling factor
    norm_factor = 1./range(curr_PPG_WK.v);
    % normalise
    wave1.v = (curr_PPG_WK.v-min(curr_PPG_WK.v)).*norm_factor;
    clear P P_out Q curr_PPG_WK norm_factor
    wave1.v = noisegen(wave1.v,SNR,seed(1));
    wave1.t = data.path_waves.aorta_finger(sim_no).onset_time(41);
    
%     wave2.v = data.path_waves.aorta_foot(sim_no).P{1, 72};  % seg. 49-ankle/anterior tibial artery
    P.v = data.path_waves.aorta_foot(sim_no).P{1,72};
    P.fs = data.path_waves.fs;
    P_out.v = data.config.p_out(sim_no,:)*ones(size(P.v));
    P_out.fs = data.path_waves.fs  ;
    Q.v = data.path_waves.aorta_foot(sim_no).U{72}.*data.path_waves.aorta_foot(sim_no).A{72};
    Q.fs = data.path_waves.fs;
    curr_PPG_WK = estimate_ppg_using_windkessel(P, Q, P_out);
    curr_PPG_WK.v = curr_PPG_WK.v(:);
    [~,rel_el] = min(curr_PPG_WK.v);
    curr_PPG_WK.v = [curr_PPG_WK.v(rel_el:end); 
    curr_PPG_WK.v(1:rel_el-1)];
    % find normalisation scaling factor
    norm_factor = 1./range(curr_PPG_WK.v);
    % normalise
    wave2.v = (curr_PPG_WK.v-min(curr_PPG_WK.v)).*norm_factor;
    clear P P_out Q curr_PPG_WK norm_factor
    wave2.v = noisegen(wave2.v,SNR,seed(2)); 
    wave2.t = data.path_waves.aorta_foot(sim_no).onset_time(72);
    
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
%%
filename = sprintf('ftPWV_n%i_%i.mat',[SNR,time]);
save(filename,'ftPWV')
%% cPP & bPP
fprintf('Add noise on PP\n')
cPP = zeros(length(data.haemods),1);
bPP = zeros(length(data.haemods),1);
for sim_no = 1:length(data.haemods)
rng(sim_no)
seed = round(rand(1,2)*10000+time*SNR);
P_bra = noisegen(data.path_waves.aorta_finger(sim_no).P{1, 24},SNR,seed(1));
P_aor = noisegen(data.path_waves.aorta_finger(sim_no).P{1, 1},SNR,seed(2));
cPP(sim_no) = max(P_aor) - min(P_aor);
bPP(sim_no) = max(P_bra) - min(P_bra);
end
filename = sprintf('cPPbPP_n%i_%i.mat',[SNR,time]);
save(filename,'cPP','bPP')
%% Pb RM aortic root
fprintf('Add noise on Pb RM\n')
wave_a.Pf_A = zeros(length(data.haemods),1);
wave_a.Pb_A = wave_a.Pf_A;
wave_a.RM = wave_a.Pf_A;
for sim_no = 1 : length(data.haemods)
    rng(sim_no)
    seed = round(rand(1,2)*10000+time*SNR);
    wave_a.p = noisegen(data.path_waves.aorta_finger(sim_no).P{1, 1},SNR,seed(1));  % seg. 1-aorta
    wave_a.v = noisegen(data.path_waves.aorta_finger(sim_no).U{1, 1},SNR,seed(2)); % seg. 1-aorta
    wave_a.rho = data.config.constants.rho(sim_no);
    wave_a.c = sqrt(sum(diff(133.322*wave_a.p).^2)./sum(diff(wave_a.v).^2))/wave_a.rho;
    wave_a.dP_f = (diff(133.322*wave_a.p) + (wave_a.rho*wave_a.c*diff(wave_a.v)))/2;
    wave_a.dP_b = (diff(133.322*wave_a.p) - (wave_a.rho*wave_a.c*diff(wave_a.v)))/2;
    wave_a.P_f = cumsum(wave_a.dP_f)/133.322;
    wave_a.P_b = cumsum(wave_a.dP_b)/133.322;
    wave_a.Pf_A(sim_no) = max(wave_a.P_f);
    wave_a.Pb_A(sim_no) = max(wave_a.P_b);
    wave_a.RM(sim_no) = wave_a.Pb_A(sim_no)/wave_a.Pf_A(sim_no);
end
Pb = wave_a.Pb_A;
RM = wave_a.RM;
filename = sprintf('PbRM_n%i_%i.mat',[SNR,time]);
save(filename,'Pb','RM')
%% AIx AP
fprintf('Add noise on AIx AP\n')
options.do_plot = 0;
options.exclude_low_quality_data = 1;
AIx = zeros(length(data.haemods),1);
AP = zeros(length(data.haemods),1);
for sim_no = 1:length(data.haemods)
    rng(sim_no)
    seed = round(rand(1,1)*10000+time*SNR);
    S.v = noisegen(data.path_waves.aorta_finger(sim_no).P{1, 1},SNR,seed);
    S.fs = 500;
    options.do_plot = 0; sig.ht = 1.75;
    [~, fid_pts, ~, ~] = PulseAnalyse10(S, options);
    p1in = fid_pts.p1in;
    p2pk = fid_pts.p2pk;
    if isnan(p1in)||isnan(p2pk)
    else
    AP(sim_no) = S.v(p2pk) - S.v(p1in);
    PP = max(S.v) - min(S.v);
    AIx(sim_no) = AP(sim_no)/PP*100;
    end
end
%% 
filename = sprintf('AIx_AP_n%i_%i.mat',[SNR,time]);
save(filename,'AP','AIx')
%% RI ppg
fprintf('Add noise on RI ppg\n')
load([PATH_ROOT,'OneDrive - King''s College London\PhD Work\code\PW_studies\pwdb_1\PWs_mat\mat\PWs_Digital.mat'])
PWs_indfinger = PWs;
RI_PPG = zeros(length(data.haemods),1);
for i = 1: length(PWs_indfinger.PPG)
    rng(sim_no)
    seed = round(rand(1,1)*10000+time*SNR);
    S.v = noisegen(PWs_indfinger.PPG{1, i},SNR,seed);
    S.fs = 500;
    options.do_plot = 0; sig.ht = 1.75;
    [~, fid_pts, ~, ~] = PulseAnalyse10(S, options);
    dia = fid_pts.dia;
    if isnan(dia)
    else
    RI_PPG(i) = S.v(dia)/max(S.v);
    end
end
%%
filename = sprintf('RI_PPG_n%i_%i.mat',[SNR,time]);
save(filename,'RI_PPG')
%% CAVI
fprintf('Add noise on CAVI\n')
curr_ind = zeros(length(data.haemods),1);
ha_pwv = curr_ind; % heart-to-ankle PWV
ha_beta = curr_ind;
a = curr_ind;
b = curr_ind;
ind_1 = [0.850 0.658 0.432];
ind_2 = [0.695 2.103 4.441];

for sim_no = 1 : length(data.haemods)
    rng(sim_no)
    seed = round(rand(1,3)*10000+time*SNR);
    P_bra = noisegen(data.path_waves.aorta_finger(sim_no).P{1, 24},SNR,seed(1));
    Ps = max(P_bra);
    Pd = min(P_bra);
    Pp = max(P_bra) - min(P_bra);
    wave1.v = noisegen(data.path_waves.aorta_foot(sim_no).P{1, 1},SNR,seed(2));   % seg. 01-aorta
    wave1.t = data.path_waves.aorta_foot(sim_no).onset_time(1);
    wave2.v = noisegen(data.path_waves.aorta_foot(sim_no).P{1, 72},SNR,seed(3));  % seg. 49-ankle/anterior tibial artery
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
ha_beta(sim_no) = log(Ps/Pd)*2*data.config.constants.rho(sim_no)*ha_pwv(sim_no)^2/(Pp*133.3224);% ln(mmHg*mmHg-1)*(kg*m-3)*(m2*s-2)/mmHg
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

curr_ind(sim_no) = a(sim_no)*2*data.config.constants.rho(sim_no)*log(Ps/Pd)*ha_pwv(sim_no)^2/(Pp*133.3224) + b(sim_no);
end
%%
CAVI = curr_ind;
filename = sprintf('CAVI_n%i_%i.mat',[SNR,time]);
save(filename,'CAVI')
%%
fprintf('Add noise on DC\n')
DC = zeros(length(data.haemods),1);
for sim_no = 1 : length(data.haemods)
    rng(sim_no)
    seed = round(rand(1,2)*10000+time*SNR);
    LA = noisegen(data.path_waves.aorta_brain(sim_no).A{1, 10},SNR,seed(1)); %carotid 10 
    PW = noisegen(data.path_waves.aorta_finger(sim_no).P{1, 24},SNR,seed(2));%brachial 24
%     LD = sqrt(LA/pi);
    LAd = min(LA);
    LAs = max(LA);
    SBP = max(PW);
    DBP = min(PW);
    PP = (SBP - DBP)*0.13332;
    DC(sim_no) = 1e3*(LAs - LAd)/(LAd*PP);
end
%%
filename = sprintf('DC_n%i_%i.mat',[SNR,time]);
save(filename,'DC')
%% Selected

age = struct2table(data.haemods).age;
clear all_var 
VA_index = table(age,aoPWV,cfPWV,baPWV,cbPWV,crPWV,ftPWV,cPP,bPP,AP,Pb,AIx,RM,RI_PPG,CAVI,DC);
var_name = VA_index.Properties.VariableNames;
var_unit = {' [years]',' [m/s]',' [m/s]',' [m/s]',' [m/s]',' [m/s]',' [m/s]',' [mmHg]',' [mmHg]',' [mmHg]',' [mmHg]',' [au]',' [au]',' [au]'...
    ,' [au]',' [10^{-3}kPa^{-1}]'};
VA_array = table2array(VA_index);
all_var.v = VA_array(data.plausibility.plausibility_log,2:end);
all_var.unit = var_unit(2:end);
all_var.name = var_name(2:end);
all_var.name{13} = 'RI_{ppg}';
all_var.age = VA_array(data.plausibility.plausibility_log,1);
filename = sprintf('noise_var_n%i_%i.mat',[SNR,time]);
save(filename,'all_var')
end
end

%%
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
E_modulus = sum((rho_weight.*((20/3)*((k(:,1).*exp(k(:,2).*radius_cm))+k(:,3))/10)),2);
clear seg_no art_len rho_weight radius radius_cm
%%
rel_rsq_s = zeros(15,4,5);
rel_rsq_p = zeros(15,4,5);
% time = 1:5;
for s=1:5
    
    snr = [30 20 15];
for no=1:3
    SNR = snr(no);
    filename = sprintf('noise_var_n%i_%i.mat',[SNR,s]);
    load(filename)
    var_arry = all_var.v;
%     eval(['all_var_n',int2str(SNR), ' = all_var;']);
    for i = 1:width(var_arry)
        refs = E_modulus;
        vals = var_arry(:,i);
        rel_els = find(~isnan(vals));
        temp = corr(refs(rel_els), vals(rel_els),'Type','Spearman');
        if isnan(temp)
            temp = 0;
        else        
        end
        rel_rsq_s(i,no+1,s) = abs(temp);
        clear temp rel_els

        rel_els = find(~isnan(vals));
        temp = corr(refs(rel_els), vals(rel_els),'Type','Pearson');
        if isnan(temp)
            temp = 0;
        else        
        end
        rel_rsq_p(i,no+1,s) = abs(temp);
        clear temp rel_els
    end
end



load('PW_indx.mat')
var_arry = table2array(PW_indx.v(data.plausibility.plausibility_log,2:end));

for i = 1:width(var_arry)
    refs = E_modulus;
    vals = var_arry(:,i);
    rel_els = find(~isnan(vals));
    temp = corr(refs(rel_els), vals(rel_els),'Type','Spearman');
    rel_rsq_s(i,1,s) = abs(temp);
    clear temp rel_els

    rel_els = find(~isnan(vals));
    temp = corr(refs(rel_els), vals(rel_els),'Type','Pearson');
    rel_rsq_p(i,1,s) = abs(temp);
    clear temp rel_els
end
end
rel_rsq_p = mean(rel_rsq_p,3);
rel_rsq_s = mean(rel_rsq_s,3);
%% Plot
curr_indx.unit = PW_indx.unit(2:end)  ;
curr_indx.name = PW_indx.name(2:end)  ;
paper_size = [900, 300];
figure('Position', [20,20,paper_size]);
var_name = string(curr_indx.name);

b = bar(rel_rsq_s,1);
set(gca,'xtick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15])
set(gca,'xticklabel',var_name)
ylim([0 1.1])
ylabel('r_s')
set(gca,'FontSize',10);

for i = 1:4
b(:,i).FaceColor = 'flat';
b(:,i).CData = 0.3*(i-1)*[1 1 1];
end
% applyhatch(gcf,'\-x/');
PATH_SAVE = [PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\figure\fig 8\'];
% saveas(gcf, [PATH_SAVE 'cc_s.png'])
%%
paper_size = [900, 300];
figure('Position', [20,20,paper_size]);
var_name = string(all_var.name);
b = bar(rel_rsq_p,1);
set(gca,'xtick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15])
set(gca,'xticklabel',var_name)
ylim([0 1.1])
ylabel('r_p')
set(gca,'XAxisLocation','top')
set(gca,'FontSize',10);

for i = 1:4
b(:,i).FaceColor = 'flat';
b(:,i).CData = 0.2*(i-1)*[1 1 1];
end

% applyhatch(gcf,'\-x/');
PATH_SAVE = [PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\figure\fig 8\'];
% saveas(gcf, [PATH_SAVE 'cc_p.png'])
%%
paper_size = [300, 600];
figure('Position', [20,20,paper_size])
fig_settings.lwidth = 1;
fig_settings.ftsize = 10;
t = tiledlayout(4,1,'TileSpacing','Compact');

nexttile
P = data.path_waves.aorta_foot(1).P{1, 1};
time = [1:length(P)]/500;
plot(time,P,'Color',0*[1 1 1],'LineWidth', fig_settings.lwidth)
ylim([72 105])
ylabel('P [mmHg]', 'FontSize', fig_settings.ftsize)
set(gca,'XTickLabel',[])
xlim([0 0.8])
title('Original')
set(gca,'FontSize', fig_settings.ftsize)
snr = [30 20 15];
for i = 1:3
    SNR = snr(i);
    nexttile
    P = noisegen(data.path_waves.aorta_foot(1).P{1, 1},SNR,2023);
    plot(time,P,'Color',0*[1 1 1],'LineWidth', fig_settings.lwidth)
    ylim([72 105])
    ylabel('P [mmHg]', 'FontSize', fig_settings.ftsize)
    if ismember(SNR,[15])
        xlabel('Time [s]', 'FontSize', fig_settings.ftsize)
    else
        set(gca,'XTickLabel',[])
    end
    tit = sprintf('SNR = %i dB',SNR);
    title(tit)
    xlim([0 0.8])
    set(gca,'FontSize', fig_settings.ftsize)
end
PATH_SAVE = [PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\figure\fig 8\'];
% saveas(gcf, [PATH_SAVE 'pressure.png'])
clearvars -except data PATH_ROOT
%%
P = data.path_waves.aorta_foot(1).P{1, 1};
Y = noisegen(P,20,0);
plot(Y);hold on
plot(P)
%%
function Y = noisegen(X,SNR,seed)
% noisegen add white Gaussian noise to a signal.
% [Y, NOISE] = NOISEGEN(X,SNR) adds white Gaussian NOISE to X. The SNR is in dB.

rng(seed);
NOISE=randn(size(X));
NOISE=NOISE-mean(NOISE);
norm = (max(X)+min(X))/2;
X = X-norm;
signal_power = 1/length(X)*sum(X.*X);
noise_variance = signal_power / (10^(SNR/10));
NOISE=sqrt(noise_variance)/std(NOISE)*NOISE;
Y=X+norm+NOISE;

end

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