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

%% Four quadrant plot
% left eta, top COE, right Im, down P
figure(1);
pBorder = 0.06;
legendVar = {};
robot_mass_void = 2183e-3; %default
leg_nbr = 6; %leg number default 
ST = 100e-3; %step length default [mm]
ERROR = 0; %basal & kinematic & basal error removal [1 0]
OPTIM = 0; %energy optimization / reduction (see for loop to set values) [1 0]
singleColor = 0; %all curves of same color [1 0]
prev_legend = 0; %use previously generated graph axis [1 0]
LDS = [1 51]%[1]; %loads indexes select %15 25 40 51

if(~prev_legend) %reset axis values
    etaMin = [];
    etaMax = [];
    cotMin = [];
    cotMax = [];
    powMin = [];
    powMax = [];
    imMin = [];
    imMax = [];
end
indicatorMask = [];
indicator = [];
indicator2 = [];

template_plot = [];
template_plot_time = [];

clc
tic
if(OPTIM), OPTIM_V = 0.7:0.01:0.98; else, OPTIM_V = 0;end %0.7:0.01:0.98
for param=OPTIM_V
for SL= [120e-3] %[70 100 120 140]*1e-3 [4 6 8] [.3 .35 .4 .45 .5 .6 .7] [.3 .5 .7]    [70 100 120 140]*1e-3 [4 6 8]
   
    switch(SL)
%           DUTY FACTOR
        case -123
            P = simu_POW_L6S110D055;
            Phip = est_POW_L6S110D055  - est_POW_L6S100_err*ERROR;
            limitsMask = isnan(est_POW_L6S110D055);
            leg_nbr = 6;
            ST  = 110e-3;
            legendVar = {legendVar{:} , 'tripod alternate gait 0.5'};
        case -1
            P = simu_POW_L6S100D05;
            Phip = est_POW_L6S100D05  - est_POW_L6S100_err*ERROR;
            limitsMask = isnan(est_POW_L6S100D05);
            leg_nbr = 6;
            ST  = 100e-3;
            legendVar = {legendVar{:} , 'tripod alternate gait 0.5'};
        case -2
            P = simu_POW_L6S100D05_wave;
            Phip = est_POW_L6S100D05_wave  - est_POW_L6S100_err*ERROR;
            limitsMask = isnan(est_POW_L6S100D05_wave);
            leg_nbr = 6;
            ST  = 100e-3;
            legendVar = {legendVar{:} , 'wave gait 0.5'};
            
        case -3
            P = simu_POW_L6S100D083_wave;
            Phip = est_POW_L6S100D083_wave  - est_POW_L6S100_err*ERROR;
            limitsMask = isnan(est_POW_L6S100D083_wave);
            leg_nbr = 6;
            ST  = 100e-3;
            legendVar = {legendVar{:} , 'wave gait 0.82'};
            
        case -4
            P = simu_POW_L4S100D06(1:20,:);
            Phip = est_POW_L4S100D06(1:20,:)  - est_POW_L4S100_err(1:20,:)*ERROR;
            limitsMask = isnan(est_POW_L4S100D06(1:20,:));
            leg_nbr = 4;
            ST  = 100e-3;
            legendVar = {legendVar{:} , 'quadripod 0.6'};
        
        case .83
            P = simu_POW_L6S100D083;
            Phip = est_POW_L6S100D083  - est_POW_L6S100_err*ERROR;
            limitsMask = isnan(est_POW_L6S100D083);
            leg_nbr = 6;
            ST  = 100e-3;
            legendVar = {legendVar{:} , '\beta = 0.83'};
            
        case .3
            P = simu_POW_L6S100D03;
            Phip = est_POW_L6S100D03  - est_POW_L6S100_err*ERROR;
            limitsMask = isnan(est_POW_L6S100D03);
            leg_nbr = 6;
            ST  = 100e-3;
            legendVar = {legendVar{:} , 'duty factor 0.3'};
        case .35
            P = simu_POW_L6S100D035;
            Phip = est_POW_L6S100D035 - est_POW_L6S100_err*ERROR;
            limitsMask = isnan(est_POW_L6S100D035);
            leg_nbr = 6;
            ST  = 100e-3;
            legendVar = {legendVar{:} , '\beta = 0.35'};
        case .4
            P = simu_POW_L6S100D04;
            Phip = est_POW_L6S100D04 - est_POW_L6S100_err*ERROR;
            limitsMask = isnan(est_POW_L6S100D04);
            leg_nbr = 6;
            ST  = 100e-3;
            legendVar = {legendVar{:} , '\beta = 0.4'};
        case .45
            P = simu_POW_L6S100D045;
            Phip = est_POW_L6S100D045 - est_POW_L6S100_err*ERROR;
            limitsMask = isnan(est_POW_L6S100D045);
            leg_nbr = 6;
            ST  = 100e-3;
            legendVar = {legendVar{:} , '\beta = 0.45'};
        case .5 
            P = simu_POW_L6S100D05;
            Phip = est_POW_L6S100D05 - est_POW_L6S100_err*ERROR;
            limitsMask = isnan(est_POW_L6S100D05);
            leg_nbr = 6;
            ST  = 100e-3;
            legendVar = {legendVar{:} , '\beta = 0.5'};
        case .6
            P = simu_POW_L6S100D06;
            Phip = est_POW_L6S100D06 - est_POW_L6S100_err*ERROR;
            limitsMask = isnan(est_POW_L6S100D06);
            leg_nbr = 6;
            ST  = 100e-3;
            legendVar = {legendVar{:} , '\beta = 0.6'};
        case .7
            P = simu_POW_L6S100D07;
            Phip = est_POW_L6S100D07 - est_POW_L6S100_err*ERROR;
            limitsMask = isnan(est_POW_L6S100D07);
            leg_nbr = 6;
            ST  = 100e-3;
            legendVar = {legendVar{:} , '\beta = 0.7'};
        
%          LEGS COUNT
        case 4
            P = simu_POW_L4S100D05;
            Phip = est_POW_L4S100D05 - est_POW_L4S100_err*ERROR;            
            P = P(1:20,:);
            Phip = Phip(1:20,:);
            limitsMask = isnan(est_POW_L4S100D05(1:20,:));
            leg_nbr = 4;
            ST  = 100e-3;
            legendVar = {legendVar{:} , 'quadripod'};
        case 6
            P = simu_POW_L6S100D05;
            Phip = est_POW_L6S100D05 - est_POW_L6S100_err*ERROR;
            limitsMask = isnan(est_POW_L6S100D05);
            leg_nbr = 6;
            ST  = 100e-3;
            legendVar = {legendVar{:} , 'hexapod'};
        case 8
            P = simu_POW_L8S100D05;
            Phip = est_POW_L8S100D05 - est_POW_L8S100_err*ERROR;
            limitsMask = isnan(est_POW_L8S100D05);
            leg_nbr = 8;
            ST  = 100e-3;
            legendVar = {legendVar{:} , 'octopod'};
            
%           STEP LENGTH            
        case 70e-3 
           P = simu_POW_L6S70D05;
           Phip = est_POW_L6S70D05 - est_POW_L6S70_err*ERROR;
           limitsMask = isnan(est_POW_L6S70D05);
            ST  = 70e-3;
            leg_nbr = 6;
            legendVar = {legendVar{:} , 'step 70mm'};
        case 100e-3
           P = simu_POW_L6S100D05;
           Phip = est_POW_L6S100D05 - est_POW_L6S100_err*ERROR;
           limitsMask = isnan(est_POW_L6S100D05);
            ST  = 100e-3;
            leg_nbr = 6;
            legendVar = {legendVar{:} , 'step 100mm'};
         case 110e-3
           P = simu_POW_L6S110D055;
           Phip = est_POW_L6S110D05 - est_POW_L6S100_err*ERROR;
           limitsMask = isnan(est_POW_L6S110D05);
           ST  = 110e-3;
           leg_nbr = 6;
           legendVar = {legendVar{:} , 'step 100mm'};
        case 120e-3
           P = simu_POW_L6S120D05;
           Phip = est_POW_L6S120D05 - est_POW_L6S120_err*ERROR;
           limitsMask = isnan(est_POW_L6S120D05);
           leg_nbr = 6;
            ST  = 120e-3;
            legendVar = {legendVar{:} , 'step 120mm'};
        case 140e-3
            P = simu_POW_L6S140D05;
            Phip = est_POW_L6S140D05 - est_POW_L6S140_err*ERROR;
            limitsMask = isnan(est_POW_L6S140D05);
            leg_nbr = 6;
            ST  = 140e-3;
            leg_nbr = 6;
            legendVar = {legendVar{:} , 'step 140mm'};
    end

PowClean = reshape(Phip(~limitsMask),[size(limitsMask,1), sum(~limitsMask(1,:))]);
fprintf("Power/load variation %f \n",mean(abs(PowClean(1,:)-PowClean(end,:))))
Im = simu_X*ST*2;
if(SL==-4 || SL==4), Im = Im(1:20,:);end
Im(limitsMask) = NaN;
Phip(limitsMask) = NaN;
P(limitsMask) = NaN;

if(OPTIM_V)
    nop=0;
    for ki=1:size(Phip,1)  
        if(min(Phip(ki,:) - min(Phip(ki,:))*param)>min(P(ki,:)))
            Phip(ki,:) = Phip(ki,:) - min(Phip(ki,:))*param; %min(Phip(ki,:))*param;

        else
          nop=1;
          break;
        end
    end
    if(nop==1)
       continue 
    end
end
%Phip = Phip - min(min(Phip))*param; %1.4*leg_nbr -  %;min(min(Phip))*0.96;
Phi_ = (Phip-P);
COE = Phi_./Im;
eth = P./Phip;

%plot a graph for each load
for i = 1:4
      h(i) = subplot(2,2,i);
end

for ld=LDS%size(simu_POW,1)  10 21
    if(numel(LDS)~=1 && ~OPTIM)
        if(ld==LDS(1))
             legendSave = legendVar{end};
             legendVar{end} = [legendSave ' ' num2str(ld)];    
        else
            legendVar = {legendVar{:} , [legendSave ' ' num2str(ld)]};
        end
    end
    
    
    % function of v (Im)
    vPhip = Phip(ld,:);
    vPhi_ = Phi_(ld,:);
    currM = (robot_mass_void+simu_Y(ld,1));
    vCOE = smooth(COE(ld,:))'/currM; %/currM
    
    vCOE(isnan(COE(ld,:))) = NaN(1,sum(isnan(COE(ld,:))));
    vP = P(ld,:);
    sIm = Im(ld,:);
    sEth = smooth(eth(ld,:))';
    sEth(isnan(eth(ld,:))) = NaN(1,sum(isnan(eth(ld,:))));
    
    template_plot = [template_plot; vCOE];
    template_plot_time = [template_plot_time; sIm];
    
    vP = vP(~isnan(vP));
    vPhip = vPhip(~isnan(vPhip));
    vPhi_ = vPhi_(~isnan(vPhi_));
    sIm = sIm(~isnan(sIm));
    sEth = sEth(~isnan(sEth));
    vCOE = vCOE(~isnan(vCOE));
    
    
    %Indicator building
    indicatorValue =min(sEth)/currM;% max(vPhi_./vPhip);%max(Phi_(ld,:)./Phip(ld,:));%sEth(1)/sEth(end); vPhi_./vPhip
    indicatorValue2 = currM;
    n1 = find(sEth==max(sEth));
    n2 = find(vP==max(vP));
    n3 = find(vCOE==min(vCOE));
    
    fprintf([num2str(n1) ' ' num2str(n2) ' ' num2str(n3) '\n']);
    
    indicator = [indicator indicatorValue];% mean(Phi_(ld,:)./Phip(ld,:))
    indicator2 = [indicator2 indicatorValue2];% mean(Phi_(ld,:)./Phip(ld,:))
    indicatorMask = [indicatorMask length(unique([n1 n2 n3]))]; 
    
%     if(sEth(1)>min(sEth) && sEth(end)>min(sEth)) %min(sEth) && sEth(end)>min(sEth) && abs(sEth(1)-sEth(end))<0.05) || (sEth(1)<max(sEth) && sEth(end)<max(sEth)) && abs(sEth(1)-sEth(end))<0.05) %sEth(1)>min(sEth) && sEth(end)>min(sEth)  abs(sEth(1)-sEth(end))<0.01
%         indicatorMask = [indicatorMask 1];
%         indicator = [indicator indicatorValue];% mean(Phi_(ld,:)./Phip(ld,:))
%     else
%         indicatorMask = [indicatorMask 0];
%         indicator = [indicator indicatorValue];
%     end
    
    subplot(221);
    hold on;
    plot(sEth,vCOE,'LineWidth',2)
    hold off;
    
    subplot(222)
    hold on;
    plot(sIm,vCOE,'LineWidth',2);
    hold off;
    
    subplot(223)
    hold on;    
    plot(sEth,vP,'LineWidth',2);    
    hold off;
    
    subplot(224)
    hold on;
    plot(sIm,vP,'LineWidth',2);
    hold off;
    
    %check axis limits
    if(~prev_legend)
        etaMin = min([etaMin sEth]);
        etaMax = max([etaMax sEth]);
        cotMin = min([vCOE cotMin]);
        cotMax = max([vCOE cotMax]);
        powMin = min([powMin vP]);
        powMax = max([powMax vP]);
        imMin = min([imMin sIm]);
        imMax = max([imMax sIm]);
    end
end

end

subplot(221)
set(h(1), 'Xdir', 'reverse')
set(h(1), 'YAxisLocation', 'Right')
set(h(1),'Position',[0+pBorder 0.5 0.5-pBorder 0.5-pBorder])
set(h(1),'xticklabel',[])
set(h(1),'yticklabel',[])
grid;
grid minor;
legend(legendVar);
legend('off')

subplot(222)
set(h(2),'Position',[0.5 0.5 0.5-pBorder 0.5-pBorder])
set(h(2),'xticklabel',[])
grid;
grid minor;
legend(legendVar,'FontSize',20);
if(singleColor || OPTIM), legend('off');end


subplot(223)
set(h(3), 'Ydir', 'reverse')
set(h(3), 'Xdir', 'reverse')
set(h(3), 'YAxisLocation', 'Right')
set(h(3),'Position',[0+pBorder 0+pBorder 0.5-pBorder 0.5-pBorder])
set(h(3),'yticklabel',[])
grid;
grid minor;
legend(legendVar);
legend('off')

subplot(224)
set(h(4), 'Ydir', 'reverse')
set(h(4),'Position',[0.5 0+pBorder 0.5-pBorder 0.5-pBorder])
grid;
grid minor;
legend(legendVar);
legend('off')

% t1 = annotation('textbox',[1-pBorder 0.5 0.001 0.001],'String','v','FontSize',12,'FontName','Calibri','EdgeColor',[1 1 1],'FontWeight','bold');
% t2 = annotation('textbox',[pBorder 0.5 0.001 0.001],'String','\eta','FontSize',12,'FontName','Calibri','EdgeColor',[1 1 1],'FontWeight','bold');
% t3 = annotation('textbox',[0.5 1-pBorder 0.001 0.001],'String','COT/M','FontSize',12,'FontName','Calibri','EdgeColor',[1 1 1],'FontWeight','bold');
% t4 = annotation('textbox',[0.5 pBorder 0.001 0.001],'String','P_{M}','FontSize',12,'FontName','Calibri','EdgeColor',[1 1 1],'FontWeight','bold');


%remake scales
h1_ticks_Y = get(h(1),'Ytick');
h2_ticks_Y = get(h(2),'Ytick');
h3_ticks_Y = get(h(3),'Ytick');
h4_ticks_Y = get(h(4),'Ytick');

h1_ticks_X = get(h(1),'Xtick');
h2_ticks_X = get(h(2),'Xtick');
h3_ticks_X = get(h(3),'Xtick');
h4_ticks_X = get(h(4),'Xtick');

h1_ticks_LY = get(h(1),'Yaxis');
h2_ticks_LY = get(h(2),'Yaxis');
h3_ticks_LY = get(h(3),'Yaxis');
h4_ticks_LY = get(h(4),'Yaxis');

h1_ticks_LX = get(h(1),'Xaxis');
h2_ticks_LX = get(h(2),'Xaxis');
h3_ticks_LX = get(h(3),'Xaxis');
h4_ticks_LX = get(h(4),'Xaxis');

if(~prev_legend)
    etaMin = etaMin *0.9;
    etaMax = etaMax * 1.1;
    imMin = imMin *0.9;
    imMax = imMax * 1.1;
    cotMin = cotMin * 0.9;
    cotMax = cotMax * 1.1;
    powMin =  powMin * 0.95;
    powMax = powMax * 1.05;
end

set(h1_ticks_LX,'Limits',[etaMin etaMax]);
set(h2_ticks_LX,'Limits',[imMin imMax]);
set(h3_ticks_LX,'Limits',[etaMin etaMax]);
set(h4_ticks_LX,'Limits',[imMin imMax]);

set(h1_ticks_LY,'Limits',[cotMin cotMax]);
set(h2_ticks_LY,'Limits',[cotMin cotMax]);
set(h3_ticks_LY,'Limits',[powMin powMax]);
set(h4_ticks_LY,'Limits',[powMin powMax]);

try
    set(h(1),'Ytick',round(cotMin:(cotMax-cotMin)/5:cotMax,0));
    set(h(2),'Ytick',round(cotMin:(cotMax-cotMin)/5:cotMax,0));
    set(h(3),'Ytick',round(powMin:(powMax-powMin)/5:powMax,0));
    set(h(4),'Ytick',round(powMin:(powMax-powMin)/5:powMax,0));
catch Me
    set(h(1),'Ytick',round(cotMin:(cotMax-cotMin)/5:cotMax,1));
    set(h(2),'Ytick',round(cotMin:(cotMax-cotMin)/5:cotMax,1));
    set(h(3),'Ytick',round(powMin:(powMax-powMin)/5:powMax,1));
    set(h(4),'Ytick',round(powMin:(powMax-powMin)/5:powMax,1));
end
    
    

set(h(1),'Xtick',round(etaMin:(etaMax-etaMin)/5:etaMax,2));
set(h(2),'Xtick',round(imMin:(imMax-imMin)/5:imMax,2));
set(h(3),'Xtick',round(etaMin:(etaMax-etaMin)/5:etaMax,2));
set(h(4),'Xtick',round(imMin:(imMax-imMin)/5:imMax,2));


set(h(1),'FontSize',14);
set(h(2),'FontSize',14);
set(h(3),'FontSize',14);
set(h(4),'FontSize',14);

for ij=1:4
    axesHandlesToAllLines = findobj(h(ij), 'Type', 'line');
    LenLeg = length(axesHandlesToAllLines);
    if(singleColor==0)
        legendColor = vz_color_generateLegend(LenLeg)./1.2;
    else
        legendColor = ones(LenLeg,3).*[99 78 31]/100;
    end
    for line=1:length(axesHandlesToAllLines)
        axesHandlesToAllLines(line).Color = legendColor(line,:);
    end
end

end
toc

%fitIndicator =  fit(indicator',indicatorMask','smoothingspline');
figure;
hold on;
%plot(fitIndicator(min(indicator):max(indicator)))
plot(indicator,indicatorMask,'or','MarkerFaceColor','red')
hold off
title('Working points = f(indicator)')
figure;
plot(indicator2,indicatorMask,'or','MarkerFaceColor','blue')
title('Working points = f(indicator2)')
%close all
figure;
plot3(indicator2,indicator,indicatorMask,'or')
title('Working points = f(indicator1,indicator2)')


indik2 = indicator2(indicatorMask==2);
indik = indicator(indicatorMask==2);

Xi = [];
Yi = [];
for i=unique(indik2)
   Yi = [Yi min(indik(indik2==i))];
   Xi = [Xi i];
end

Xi
Yi