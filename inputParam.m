%% Theoretical PWV - baseline 
sim_no = 1;
A = [0.000, 0.000, 0.000, 0.000, 0.044];
B = [0.83, 0.99, 1.05, 1.18, 0.85]*1e-3;
C = [5.55, 5.69, 5.91, 6.17, 5.73];
sbp = [115, 125, 135, 150, 170];
dbp = [77.5, 82.5, 87.5, 95, 105];
init_mbp_vals = dbp+(0.4*(sbp-dbp));  % calculate approximate MBP corresponding to these SBP and DBP values using 0.4 constant from article; similar to 0.412 constant from: "Formula and nomogram for the sphygmomanometric calculation of the mean arterial pressure", http://heart.bmj.com/content/heartjnl/84/1/64.full.pdf

% extend range of MBP values by extrapolating at:
% - lower end
no_to_interpolate = 5;
interp_els = (1-no_to_interpolate):0;
new_sbp_vals = interp1(1:3, sbp(1:3), interp_els, 'linear', 'extrap');
new_dbp_vals = interp1(1:3, dbp(1:3), interp_els, 'linear', 'extrap');
new_A_vals = interp1(1:3, A(1:3), interp_els, 'linear', 'extrap');
new_B_vals = interp1(1:3, B(1:3), interp_els, 'linear', 'extrap');
new_C_vals = interp1(1:3, C(1:3), interp_els, 'linear', 'extrap');
% store these new values
sbp = [new_sbp_vals, sbp];
dbp = [new_dbp_vals, dbp];
A = [new_A_vals, A];
B = [new_B_vals, B];
C = [new_C_vals, C];

% - upper end
no_to_interpolate = 3;
interp_els = 2+(1:no_to_interpolate);
new_sbp_vals = interp1(1:2, sbp([end-1, end]), interp_els, 'linear', 'extrap');
new_dbp_vals = interp1(1:2, dbp([end-1, end]), interp_els, 'linear', 'extrap');
new_A_vals = interp1(1:2, A([end-1, end]), interp_els, 'linear', 'extrap');
new_B_vals = interp1(1:2, B([end-1, end]), interp_els, 'linear', 'extrap');
new_C_vals = interp1(1:2, C([end-1, end]), interp_els, 'linear', 'extrap');
% store these new values
sbp = [sbp, new_sbp_vals];
dbp = [dbp, new_dbp_vals];
A = [A, new_A_vals];
B = [B, new_B_vals];
C = [C, new_C_vals];

% calculate mbp
mbps.vals = dbp+(0.4*(sbp-dbp));  % approximation as above: 0.4 constant from article; 0.412 constant from: http://heart.bmj.com/content/heartjnl/84/1/64.full.pdf
mbps.inds = 1 : length(mbps.vals);

% data from Mattace-Raso2010 table 5 (median and 10 and 90 pc, therefore these represent +/- 40% in each direction)
no_sds = range(norminv([0.5, 0.9])); % in this case, for 80% confidence interval
%%
sd_percentile.lower_v = (1/no_sds)*100* [(6.0-5.2)/6.0, (6.4-5.7)/6.4, (6.7-5.8)/6.7, (7.2-5.7)/7.2, (7.6-5.9)/7.6; ...
    (6.5-5.4)/6.5, (6.7-5.3)/6.7, (7.0-5.5)/7.0, (7.2-5.5)/7.2, (7.6-5.8)/7.6; ...
    (6.8-5.8)/6.8, (7.4-6.2)/7.4, (7.7-6.5)/7.7, (8.1-6.8)/8.1, (9.2-7.1)/9.2; ...
    (7.5-6.2)/7.5, (8.1-6.7)/8.1, (8.4-7.0)/8.4, (9.2-7.2)/9.2, (9.7-7.4)/9.7; ...
    (8.7-7.0)/8.7, (9.3-7.6)/9.3, (9.8-7.9)/9.8, (10.7-8.4)/10.7, (12.0-8.5)/12.0; ...
    (10.1-7.6)/10.1, (11.1-8.6)/11.1, (11.2-8.6)/11.2, (12.7-9.3)/12.7, (13.5-10.3)/13.5 ...
    ];
sd_percentile.upper_v = (1/no_sds)*100* [(7.0-6.0)/6.0, (7.5-6.4)/6.4, (7.9-6.7)/6.7, (9.3-7.2)/7.2, (9.9-7.6)/7.6; ...
    (7.9-6.5)/6.5, (8.2-6.7)/6.7, (8.8-7.0)/7.0, (9.3-7.2)/7.2, (11.2-7.6)/7.6; ...
    (8.5-6.8)/6.8, (9.0-7.4)/7.4, (9.5-7.7)/7.7, (10.8-8.1)/8.1, (13.2-9.2)/9.2; ...
    (9.2-7.5)/7.5, (10.4-8.1)/8.1, (11.3-8.4)/8.4, (12.5-9.2)/9.2, (14.9-9.7)/9.7; ...
    (11.4-8.7)/8.7, (12.2-9.3)/9.3, (13.2-9.8)/9.8, (14.1-10.7)/10.7, (16.5-12.0)/12.0; ...
    (13.8-10.1)/10.1, (15.5-11.1)/11.1, (15.8-11.2)/11.2, (16.7-12.7)/12.7, (18.2-13.5)/13.5 ...
    ];
%%
% add on extrapolated values (assuming a constant percentage value for the SD outside of the specified data range)
sd_percentile.mbp = [0, init_mbp_vals, 500];
sd_percentile.lower_v = [sd_percentile.lower_v(:,1), sd_percentile.lower_v, sd_percentile.lower_v(:,end)];
sd_percentile.upper_v = [sd_percentile.upper_v(:,1), sd_percentile.upper_v, sd_percentile.upper_v(:,end)];

% find reference PWV values for different ages and different MBPs
age = 20:80;
for n = 1 : length(A)
    pwv_age(:,n) = A(n)*age + B(n)*(age.^2) + C(n);
end
%%
sims_co = data.config.hr.*data.config.sv/1000; % in l/min
mean_flow = sims_co/(1000*60); % in m3/sec
network_r = (133.33*data.config.constants.p_drop)./mean_flow;
sims_mbp = ((mean_flow.*(data.config.pvr+network_r)) + data.config.p_out)/133.33; 
%% age-specific
curr_age = data.config.age(1);
[~, rel_age_el] = min(abs(age-curr_age));
rel_pwv_age_vals = pwv_age(rel_age_el,:);
%% mbp specific

mbps.vals = dbp+(0.4*(sbp-dbp));  % approximation as above: 0.4 constant from article; 0.412 constant from: http://heart.bmj.com/content/heartjnl/84/1/64.full.pdf
mbps.inds = 1 : length(mbps.vals);
 rel_mbp_ind = interp1(mbps.vals, mbps.inds, sims_mbp(1), 'linear', 'extrap');
 %%
if rel_mbp_ind > length(mbps.vals)
    rel_mbp_ind = length(mbps.vals);
end
if rel_mbp_ind < 1
    rel_mbp_ind = 1;
end