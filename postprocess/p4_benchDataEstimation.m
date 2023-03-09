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

%% Estimate power at fixed beta (variable load and speed)

%SET GAIT PARAMETERS
step_length = 100e-3; %Not used, just for information
duty_factor = 0.5;
GaitList;gaitSelect = gait_tripod;COT_INVERT_GAIT = 0; %0 : normal %gait_tripod
legs_nbr = length(gaitSelect);
robot_mass =legs_nbr*data.mass  + 0.38; 

%MODEL LIMITS
COT_DISPLAY_LIMIT_UP = 200;
COT_DISPLAY_LIMIT_DOWN = 0;
MODEL_ACCEPTANCE = 0; %[s] : how far we can go outside experimental data
COT_RESOLUTION = 500;

%Check model boundaries
boundStmin = min(expStancePeriod_info)-MODEL_ACCEPTANCE;
boundStmax = max(expStancePeriod_info(~isinf(expStancePeriod_info)))+MODEL_ACCEPTANCE;
boundSwmin = min(expSwingPeriod_info(expSwingPeriod_info>0))-MODEL_ACCEPTANCE;
boundSwmax = max(expSwingPeriod_info) + MODEL_ACCEPTANCE;

boundLoadmin = 0;
boundLoadmax = max(expLoad_info)*legs_nbr;

boundSpeedmin = 1/(boundStmax+boundSwmax); %2*step_length/
boundSpeedmax = 1/(boundStmin+boundSwmin); %2*step_length/

robot_load = boundLoadmin:(boundLoadmax-boundLoadmin)/COT_RESOLUTION:boundLoadmax;  
robot_speed = boundSpeedmin:(boundSpeedmax-boundSpeedmin)/COT_RESOLUTION:boundSpeedmax; %m/s

[X,Y] = meshgrid(robot_speed,robot_load); %X freqs Y loads
Z_POW = NaN(size(X));

%Define gait
stepPattern = GaitPlotter(duty_factor,gaitSelect,'Resolution',COT_RESOLUTION,'T',1);
if(COT_INVERT_GAIT)                                           
  stepPattern = stepPattern(:,:,1); 
else
  stepPattern = ~stepPattern(:,:,1);
end
%Number of legs touching the ground
%study the case of a signle leg
loadDivPattern_st = sum(stepPattern).*stepPattern(1,:); %number of stance legs
loadDivPattern_sw = sum(~stepPattern).*(~stepPattern(1,:)); %number of swing legs

%Calculating each different type of power value
loadDiv_st = unique(loadDivPattern_st(loadDivPattern_st~=0));
loadDiv_sw = unique(loadDivPattern_sw(loadDivPattern_sw~=0));
loadDiv_st_count = numel(loadDiv_st);
loadDiv_sw_count = numel(loadDiv_sw);


pBar = progressBar([length(robot_load) length(robot_speed)], 5);

for ldi=1:numel(robot_load)
   pBar.update
   ld = robot_load(ldi);
   total_mass = robot_mass + ld;
    for spi=1:numel(robot_speed)
       pBar.update
       sp = robot_speed(spi);
        
       step_period = 1/sp; %2*step_length/
       stance_period = duty_factor*step_period;
       swing_period = (1-duty_factor)*step_period;
       
       Test_stance = stance_period > boundStmin && stance_period < boundStmax;
       Test_swing = swing_period > boundSwmin && swing_period < boundSwmax;

          if(Test_stance && Test_swing) %Checking limits              
              powerValues = NaN(1,loadDiv_st_count+loadDiv_sw_count);
              for i=1:loadDiv_st_count
                 currMass = (total_mass-data.mass*(loadDiv_st(i)))/loadDiv_st(i);
                 if(currMass<max(expLoad_info))
                    powerValues(i) = sum(loadDivPattern_st == loadDiv_st(i))*fitFct_st(1/stance_period,currMass);
                 end
              end
              for i=1:loadDiv_sw_count
                  if(legs_nbr==loadDiv_sw(i))
                      currMass = 0;
                  else
                      currMass = (total_mass-data.mass*(legs_nbr-loadDiv_sw(i)))/(legs_nbr-loadDiv_sw(i));
                  end
                  
                 if(currMass<max(expLoad_info))
                    powerValues(loadDiv_st_count+i) = sum(loadDivPattern_sw == loadDiv_sw(i))*fitFct_sw(1/swing_period,currMass); 
                 end
              end
              
              if(sum(isnan(powerValues))==0)
                powerValues_bk = powerValues;
                
                  leg_power = sum(powerValues)/COT_RESOLUTION;          %check this
                  robot_power = legs_nbr*leg_power;
                  if(robot_power<20)
                     disp(robot_power); 
                  end
                  %CoT = robot_power/(total_mass*GRAVITY*sp);
                  %if(CoT>COT_DISPLAY_LIMIT_DOWN && CoT<COT_DISPLAY_LIMIT_UP)
                  Z_POW(ldi,spi) = robot_power;
                  %end
              
              end
              
          end

    end
end

Z_POW = round(Z_POW,2);

figure;
surf(X,Y,Z_POW,'EdgeColor','none','FaceColor','interp')
xlabel('X - Speed [m/s]')
ylabel('Y - Load [kg]')
zlabel('Z - Power [W]') %Cost of transport
cbar = colorbar;
ylabel(cbar, 'Estimated robot power','FontSize',12) %'Estimated robot CoT'
grid on
view(0,90)


%% Make compatible with simulation values
% Limit by minimum period of each phase
%Simplify known power values
est_POW = NaN(size(simu_POW));
est_POW_mask = ones(size(simu_POW)); %Out of speed values

dx = abs(X(1,1)-X(1,2)); %estimation resolution
dy = abs(Y(1,1)-Y(2,1));

for xi=1:size(simu_X,1)
    for yi=1:size(simu_X,2)
        sx = simu_X(xi,yi); %freqs
        sy = simu_Y(xi,yi); %%loads
        
        stLim = 1/sx * duty_factor;
        swLim = 1/sx * (1-duty_factor);
        if(stLim < boundStmin || swLim < boundSwmin)
            est_POW_mask(xi,yi) = 0;
        end

        x = vz_array_findNearest(X(1,:),sx); %spd
        y = vz_array_findNearest(Y(:,1),sy); %ld
        
        if(abs(X(1,x)-sx)<=2*dx)
            if(abs(Y(y,1)-sy)<=2*dy)
                if(~isnan(Z_POW(y,x)))
                    est_POW(xi,yi) = Z_POW(y,x);
                    %fprintf('x=%f y=%f \n',sx,sy)
                end
            end
        end 
    end
end

%Fill missing by polynomial interpolation

for sColumn=1:size(est_POW,2) %vertical fit
    sYColumn = (~isnan(est_POW(:,sColumn)));
    if(sum(sYColumn)>=3 && sum(~sYColumn)~=0) %minimal number of points to do vertical fit
        sFit = fit(simu_Y(sYColumn,sColumn),est_POW(sYColumn,sColumn),'poly2');
        est_POW(~sYColumn,sColumn) = sFit(simu_Y(~sYColumn,sColumn));
    else
        continue
    end
    %figure;hold on;plot(simu_Y(sYColumn,sColumn),est_POW(sYColumn,sColumn));plot(sFit);hold off
end

for sLine=1:size(est_POW,1) %horizontal fit
    sYLine = (~isnan(est_POW(sLine,:)));
    if(sum(~sYLine)>0 && sum(sYLine)>=3) %if there is missing points
        sFit = fit(simu_X(sLine,sYLine)',est_POW(sLine,sYLine)','poly2');
        est_POW(sLine,~sYLine) = sFit(simu_X(sLine,~sYLine));
    else
        continue
    end
    %figure;hold on;plot(simu_Y(sYColumn,sColumn),est_POW(sYColumn,sColumn));plot(sFit);hold off
end

%Remove out of speed values
for xi=1:size(simu_POW,1)
    for yi=1:size(simu_POW,2)
        if(~est_POW_mask(xi,yi))
           est_POW(xi,yi) = NaN;
        end
    end
end

figure;
surf(simu_X,simu_Y,est_POW,'EdgeColor','none','FaceColor','interp')
xlabel('X - Speed [m/s]')
ylabel('Y - Load [kg]')
zlabel('Z - Power [W]') %Cost of transport
cbar = colorbar;
ylabel(cbar, 'Estimated robot power','FontSize',12) %'Estimated robot CoT'
grid on
view(0,90)

%est_POW_L6S100D07 = round(est_POW,2);
%est_POW_L4S100D05 = est_POW_L4S100D05 .*(est_POW_L4S100D05>0)

%% Extrapolate for various steps
%Estimate k values for load frequency

fitStep1 = fitFct_steps1(simu_X,simu_Y);
fitStep2 = fitFct_steps2(simu_X,simu_Y);
fitStep3 = fitFct_steps3(simu_X,simu_Y);

targetStep = 110; %step /2 in mm
est_POW_adj = est_POW.*(fitStep1*targetStep^2 + fitStep2*targetStep + fitStep3);

%Plot comparison

figure;
hold on;
surf(simu_X,simu_Y,est_POW,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.5)
surf(simu_X,simu_Y,est_POW_adj,'EdgeColor','none','FaceColor','interp')
hold on;
xlabel('X - Speed [m/s]')
ylabel('Y - Load [kg]')
zlabel('Z - Power [W]') %Cost of transport
cbar = colorbar;
ylabel(cbar, 'Estimated robot power','FontSize',12) %'Estimated robot CoT'
grid on
view(0,90)

%est_POW_L6S110D05 = round(est_POW_adj,2);

