% clear all;clc
% PATH_ROOT = 'C:\Users\43897\';
PATH_ROOT = 'C:\Users\Jingyuan Hong\';
% PATH_ROOT = 'C:\Users\jh22\';
% load([PATH_ROOT,'Downloads\pwdb_data_w_aorta_finger_path'])
load([PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\code\data_paper\aoPWV.mat'])
load([PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\code\data_paper\cfPWV.mat'])
load([PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\code\data_paper\baPWV.mat'])
load([PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\code\data_paper\cbPWV.mat'])
load([PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\code\data_paper\crPWV.mat'])
load([PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\code\data_paper\ftPWV.mat'])
load([PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\code\data_paper\PP.mat'])
load([PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\code\data_paper\AIx_AP.mat'])
load([PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\code\data_paper\Pb_RM.mat'])
load([PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\code\data_paper\RI_PPG.mat'])
load([PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\code\data_paper\CAVI.mat'])
load([PATH_ROOT,'OneDrive - King''s College London\PhD Work\2022\Comprehensive Vascular Ageing Analysis\code\data_paper\DC.mat'])
%%
clear PW_indx
age = struct2table(data.haemods).age;
PW_indx.v = table(age,aoPWV,cfPWV,baPWV,cbPWV,crPWV,ftPWV,cPP,bPP,AP,Pb ...
                ,AIx,RM,RI_PPG,CAVI,DC);
PW_indx.name = PW_indx.v.Properties.VariableNames;
PW_indx.unit = {' [years]',' [m/s]',' [m/s]',' [m/s]',' [m/s]',' [m/s]',' [m/s]',' [mmHg]',' [mmHg]',' [mmHg]',' [mmHg]',' [au]',' [au]',' [au]'...
    ,' [au]',' [10^{-3}kPa^{-1}]'};
PW_indx.name{14}=  'RI_{ppg}';
save PW_indx.mat PW_indx
clearvars -except data PATH_ROOT
%% data for Fig 4



