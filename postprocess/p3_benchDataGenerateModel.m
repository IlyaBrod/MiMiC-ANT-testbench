%% MiMiC-Ant bench programm, controls a robotic leg analyse test bench.
%     Copyright (C) 2022-2023  Ilya Brodoline - ISM AMU, Marseille (France).
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
% 
% contact: ilya.brodoline@univ-amu.fr

%% Extract data from bench recordings

SETTINGS

expNbr = numel(data.expList);   %number of recording sessions
expStancePower_info = zeros(expNbr,1);  %store power values
expSwingPower_info = zeros(expNbr,1);   %store power values
expSwingPeriod_info = zeros(expNbr,1);  %store period timings
expStancePeriod_info = zeros(expNbr,1); %store period timings
expSpeed_info = zeros(expNbr,1); %store walk speed values [%]
expLoad_info = zeros(expNbr,1); %store transported load values

fprintf("Processing \n");
pBar = progressBar(expNbr,10);
for idExp=1:expNbr
    pBar.update;
    
    currExp = data.expList(idExp); 
    expSpeed_info(idExp) = data.speedList(idExp);
    expLoad_info(idExp) = data.loadList(idExp);

    %Import data
    lMotion = currExp.legMotion;
    lvoltage = currExp.voltage;
    if(length(lvoltage)>1), lvoltage=lvoltage(:,2);end
    lcurrent = currExp.current;
    
    lForce = currExp.force;
    [~, lForce]=resampleData(lForce,lcurrent);
    lForce = abs(lForce);
    
    loc_z = lMotion(:,4);
    %loc_y = currExp.legMotion(:,3);
    loc_time = lMotion(:,1);
    power = lvoltage.*lcurrent(:,2);
    power_time = lcurrent(:,1);
    
    clear lvoltage lcurrent lMotion
    
    %case of static leg
    if(data.speedList(idExp)==0)
        expStancePower_info(idExp) = mean(power);
        expSwingPower_info(idExp) = NaN;
        expSwingPeriod_info(idExp) = 0;
        expStancePeriod_info(idExp) = Inf;
        continue;
    end
    
    %Resample data
    [~, loc_z]=resampleData([loc_time loc_z],[power_time power]);
    %[~, loc_y]=resampleData([loc_time loc_y],[power_time power]);
    %[~, loc_z]=resampleData([loc_time loc_z],[power_time power]);
    tech = 1/getSampleRate(power_time);
    
    %Apply delay & crop data
    delay = synchData(loc_z,-detrend(lForce));
    nZ = length(loc_z);
    nP = length(power);
    if(abs(delay)<nZ/3 || abs(delay)<nP/3)
        if(delay>=0)
           %loc_time = loc_time(1:nZ-delay);
           loc_z = loc_z(1+delay:end);
           %power_time = power_time(1:nZ-delay);
           power = power(1:nZ-delay);
           lForce = lForce(1:nZ-delay);
           uTime = power_time(1:nZ-delay);
        else
           %power_time = power_time(1:nP+delay);
           uTime = power_time(1:nP+delay);
           power = power(1-delay:end);
           lForce = lForce(1-delay:end);
           %loc_time = loc_time(1:nP+delay);
           loc_z = loc_z(1:nP+delay);
        end
    else
        if(nZ>nP)
           loc_z = loc_z(1:length(nP));
           uTime = power_time;
        else
           power = power(1:length(nZ));
           lForce = lForce(1:length(nZ));
           uTime = power_time(1:length(nZ));
        end
    end
    
    uTime = (1:length(loc_z))*tech;
    
    %DEBUG
    %figure; hold on;plot(uTime,lForce);plot(uTime,loc_z);hold off;
    %figure; hold on;plot(uTime,lForce);plot(uTime,power);hold off;
    
    %Calculate power
        error = COT_GROUND_ERROR;

        %detect ground contact
        %rloc_z = rescale(loc_z);
        [modx,mody] = modeMap(loc_z,error);
        [~,locs] = findpeaks(mody,'MinPeakProminence',max(mody)/2);

        %split stance / swing timing
        groundZ = min(modx(locs)); %Can remplace by known value from calibration
        stance_idx = loc_z< groundZ + error & loc_z > groundZ - error;
        stance_time = uTime(stance_idx); %same as loc_time resampled
        %rstance_time = (0:length(stance_time)-1)'*tech; %resample
        swing_idx = ~stance_idx;
        swing_time = uTime(swing_idx);
        %rswing_time = (0:length(swing_time)-1)'*tech;

        %power calculation
        stance_power = power(stance_idx);
        swing_power = power(swing_idx);
    
        plotPower_swing = NaN(size(power));
        plotPower_swing(swing_idx) = power(swing_idx);
        plotPower_stance = NaN(size(power));
        plotPower_stance(stance_idx) = power(stance_idx);
        % DEBUG PLOTS
        %figure; hold on; plot(uTime,power); plot(stance_time,stance_power); hold off
        %figure;yyaxis left; hold on; plot(uTime,plotPower_swing); plot(uTime,plotPower_stance);yyaxis right; plot(uTime,loc_z); hold off
        %legend('swing','stance')
        %figure; hold on; plot(uTime,loc_z); plot(uTime(stance_idx),loc_z(stance_idx),'or','MarkerFaceColor','red'); hold off
        
        
        %Calculate phases timings
        swing_diff = round(diff(swing_time),3);
        %figure;hold on;plot(swing_time,loc_z(swing_idx));plot(swing_time(1:end-1),swing_diff,'or','MarkerFaceColor','red');hold off;
        [xmode_sw, ymode_sw] = modeMap(swing_diff(swing_diff>max(swing_diff)/5),10^(-COT_PERIODS_DECIMALS));
        [~,sw_locs] = findpeaks(swing_diff,'MinPeakProminence',min(xmode_sw(ymode_sw==max(ymode_sw))*0.90));   
        swing_periods = medfilt1(swing_time(sw_locs+1)-swing_time(sw_locs),COT_GROUND_MEDIAN,'truncate');
        %swing_periods = swing_periods(COT_GROUND_MEDIAN:end);
        
        stance_diff = round(diff(stance_time),3);
        %figure;hold on;plot(stance_time,loc_z(stance_idx));plot(stance_time(1:end-1),stance_diff,'or','MarkerFaceColor','red');hold off;
        [xmode_st, ymode_st] = modeMap(stance_diff(stance_diff>max(stance_diff)/5),10^(-COT_PERIODS_DECIMALS));
        [~,st_locs] = findpeaks(stance_diff,'MinPeakProminence',min(xmode_st(ymode_st==max(ymode_st))*0.90));
        stance_periods = medfilt1(stance_time(st_locs+1)-stance_time(st_locs),COT_GROUND_MEDIAN,'truncate');
        %stance_periods = stance_periods(COT_GROUND_MEDIAN:end);
        
        % Split power into steps
        steps_number = min(numel(sw_locs),numel(st_locs));
        swingPowers_temp = zeros(1,steps_number-2);
        stancePowers_temp = zeros(1,steps_number-2);
        
        for step=2:steps_number-1
                swing_power_seg = swing_power(sw_locs(step-1):sw_locs(step));
                stance_power_seg = stance_power(st_locs(step-1):st_locs(step));
                
                    %Estimate mean value
                    swingPowers_temp(step-1) = mean(swing_power_seg);
                    stancePowers_temp(step-1) = mean(stance_power_seg);
        end

        %Save data : estimate mean stance & mean swing power
        expStancePower_info(idExp) =  mean(stancePowers_temp);
        expSwingPower_info(idExp) = mean(swingPowers_temp);
        
        expSwingPeriod_info(idExp) = round(mean(swing_periods),COT_PERIODS_DECIMALS);
        expStancePeriod_info(idExp) = round(mean(stance_periods),COT_PERIODS_DECIMALS);

    clear power power_time loc_z loc_time stance_power swing_power lMotion lvoltage lcurrent
end

fprintf("Done\n")


%% Format data
xSpeed = unique(data.speedList); %for plot purcent of speed  
xTime_st = unique(expStancePeriod_info); 
xTime_sw = unique(expSwingPeriod_info);
yLoad = unique(data.loadList);
xSpeed_count = numel(xSpeed);
xTime_st_count = numel(xTime_st);
xTime_sw_count = numel(xTime_sw);
yLoad_count = numel(yLoad);

if(containsArr(xSpeed,0)), nullSpeed = 1; else, nullSpeed=0;end 

%% Fit Stance
[X_st,Y_st] = meshgrid(1./xTime_st,yLoad);
Z_st = NaN(yLoad_count,xTime_st_count);

for x=1:xTime_st_count
    [~,idx] = containsArr(expStancePeriod_info,xTime_st(x));
    powerList = expStancePower_info(idx{1});
    for y=1:yLoad_count 
        [~,idxL] = containsArr(expLoad_info(idx{1}),yLoad(y));
        if(~isempty(idxL{1}))
            Z_st(y,x) = mean(powerList(idxL{1}));
        end
    end
end
Z_st = fillmissing(Z_st,'linear',1);
Z_st = fillmissing(Z_st,'linear',2);
xstData = reshape(X_st,[size(X_st,1)*size(X_st,2),1]);
ystData = reshape(Y_st,[size(Y_st,1)*size(Y_st,2),1]);
zstData = reshape(Z_st,[size(Z_st,1)*size(Z_st,2),1]);
ft = fittype( 'loess' );
opts = fitoptions( 'Method', 'LowessFit' );
opts.Normalize = 'on';
opts.Span = 0.8;
[fitFct_st, gof_st] = fit( [xstData, ystData], zstData, ft ,opts);
plot(fitFct_st, [xstData, ystData], zstData);

%% Fit Swing power

[X_sw,Y_sw] = meshgrid(1./xTime_sw,yLoad);
Z_sw = NaN(yLoad_count,xTime_sw_count);

for x=1:xTime_sw_count
    [~,idx] = containsArr(expSwingPeriod_info,xTime_sw(x));
    powerList = expSwingPower_info(idx{1});
    for y=1:yLoad_count 
        [~,idxL] = containsArr(expLoad_info(idx{1}),yLoad(y));
        if(~isempty(idxL{1}))
            Z_sw(y,x) = mean(powerList(idxL{1}));
        end
    end
end
if(nullSpeed) 
    Z_sw = Z_sw(:,2:end);
    X_sw = X_sw(:,2:end);
    Y_sw=Y_sw(:,2:end);
end

Z_sw = fillmissing(Z_sw,'linear',1);
Z_sw = fillmissing(Z_sw,'linear',2);
xswData = reshape(X_sw,[numel(X_sw),1]);
yswData = reshape(Y_sw,[numel(Y_sw),1]);
zswData = reshape(Z_sw,[numel(Z_sw),1]);
ft = fittype( 'loess' );
opts = fitoptions( 'Method', 'LowessFit' );
opts.Normalize = 'on';
opts.Span = 0.9; %0.6
[fitFct_sw, gof_sw] = fit( [xswData, yswData], zswData, ft ,opts);
plot(fitFct_sw, [xswData, yswData], zswData);

%% Plot data


subplot(121)                    %STANCE 3D (2)
ymax = max(ystData); ymin = min(ystData);
xmax = max(xstData); xmin = min(xstData(xstData~=0));
x_plt = xmin:(xmax-xmin)/COT_PLOT_RESOLUTION:xmax;
y_plt = min(ystData):(ymax-ymin)/COT_PLOT_RESOLUTION:ymax;
x_plt = x_plt(1:min(length(x_plt),length(y_plt))); %Crop by security
y_plt = y_plt(1:min(length(x_plt),length(y_plt)));
[X,Y] = meshgrid(x_plt,y_plt);
Z = round(fitFct_st(X,Y),1);
[X,Y] = meshgrid(1./x_plt,y_plt);

surf(X,Y,Z,'EdgeColor','none','FaceColor','interp')
xlabel('Stance period T_{stance} [s]')
ylabel('Load [kg]')
zlabel('Average stance power [W]')
grid on
%set(gca,'xscale','log')
ylim([ymin ymax])
cbar = colorbar;
ylabel(cbar, 'Average stance power estimation [W]','FontSize',12)
%title(["Stance power estimation",strcat('R^2= ',sprintf('%.3f',gof.rsquare))])
view(0,90)
set(gca,'ColorMap', cmap);

subplot(122)                    %SWING 3D (2)
xlabel('Load [kg]')
ymax = max(yswData); ymin = min(yswData);
xmax = max(xswData); xmin = min(xswData);
x_plt = xmin:(xmax-xmin)/COT_PLOT_RESOLUTION:xmax;
y_plt = min(yswData):(ymax-ymin)/COT_PLOT_RESOLUTION:ymax;
[X_sw,Y_sw] = meshgrid(x_plt,y_plt);
Z = round(fitFct_sw(X_sw,Y_sw),1);
[X,Y] = meshgrid(1./x_plt,y_plt);

surf(X,Y,Z,'EdgeColor','none','FaceColor','interp')
xlabel('Swing period T_{swing} [s]')
ylabel('Load [kg]')
grid on
%set(gca,'xscale','log')
ylim([ymin ymax])
cbar = colorbar;
ylabel(cbar, 'Average swing power estimation [W]','FontSize',12)
%title(["Stance power estimation",strcat('R^2= ',sprintf('%.3f',gof.rsquare))])
view(0,90)

%set(gca,'ColorMap', cmap); import from vz_tools custom github toolbox