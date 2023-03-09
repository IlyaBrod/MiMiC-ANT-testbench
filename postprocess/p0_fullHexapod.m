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

%% Experiment parameters
path = './data'; %path to mocap data

voltage = 11.6;
robot_mass_fr = 2183;
g = 9.80665;%m/s
freqList=[0 0.5 1 1.5 2 2.5 3];
loadList=[0 500 1000 2000];
stepList = [20 35 50 60]; % = step length /2 (defined by robot control code)
Dtip = 51*1e-3; %Distance mocap marker to tip (need to add also z offset, but not obligatory)

%% Import data from experiment files

load([path 'zero.mat'])
zeroCurr = mean(current);

%load manually defined indexes
load('Zhand.mat')
Zc1b = Zc1;
Zc2b = Zc2;
load('Zhand_more.mat')
Zc1b(:,:,3:4) = Zc1(:,:,3:4);
Zc2b(:,:,3:4) = Zc2(:,:,3:4);


[X,Y,Z] = meshgrid(freqList,loadList,stepList);
%Zf = Z; %frequencies storage
%Zd = Z; %distance of step
%Zs1 = Z; %step start
%Zc1 = Z; %current start
%Zc2 = Z; %current stop
Zspd = Z;
Zpow = Z;

pBar = progressBar([numel(stepList) numel(freqList) numel(loadList)],5);
for step = 1:numel(stepList)
    for freq = 1:numel(freqList)
        for ld = 1:numel(loadList)
            pBar.update;
            %fprintf("Progress : Step %.2f FREQ %.2f LD %.2f \n",stepList(step),freqList(freq),loadList(ld));
            
            %Load energetic values
                S = load(['LD_' num2str(loadList(ld)) '_FREQ_' num2str(freqList(freq)) '_DX_' num2str(stepList(step)) '.mat']);
                power_time = S.time;
                power_curr = S.current-zeroCurr;
                power_fe = getSampleRate(power_time);
                power = power_curr*voltage;
                
                %check static robot
                if(freqList(freq)==0)
                    Z(ld,freq,step) = mean(power);
                   continue
                end
            
            %Import trajectories
                [time_body, X_body, Y_body, Z_body,Xr_body,Yr_body,Zr_body] = openPos(loadList(ld), freqList(freq),stepList(step),'antbot_body', path);
                [time_leg, X_leg, Y_leg, Z_leg,Xr_leg,Yr_leg,Zr_leg] = openPos(loadList(ld), freqList(freq),stepList(step),'antbot_leg', path);
                fe = getSampleRate(time_body);
                
                %Rescale
                [X_body,X_leg] = setSameSize(X_body,X_leg);
                [Y_body,Y_leg] = setSameSize(Y_body,Y_leg);
                [Z_body,Z_leg] = setSameSize(Z_body,Z_leg);
                [time_body,time_leg] = setSameSize(time_body,time_leg);
                
            %Filter trajectories
%                 X_body = X_body(round(numel(X_body)*0.50):end);
%                 Y_body = Y_body(round(numel(Y_body)*0.50):end);
%                 Z_body = Z_body(round(numel(Z_body)*0.50):end);

%                 X_leg = X_leg(round(numel(X_leg)*0.50):end);
%                 Y_leg = Y_leg(round(numel(Y_leg)*0.50):end);
%                 Z_leg = Z_leg(round(numel(Z_leg)*0.50):end);

                X_body = medfilt1(X_body,round(numel(X_body)*0.05),'omitnan','truncate');
                Y_body = medfilt1(Y_body,round(numel(Y_body)*0.05),'omitnan','truncate');
                Z_body = medfilt1(Z_body,round(numel(Z_body)*0.05),'omitnan','truncate');
                X_leg = medfilt1(X_leg,round(numel(X_leg)*0.05),'omitnan','truncate');
                Y_leg = medfilt1(Y_leg,round(numel(Y_leg)*0.05),'omitnan','truncate');
                Z_leg = medfilt1(Z_leg,round(numel(Z_leg)*0.05),'omitnan','truncate');
                
                %Find start and end of movement
                d_body = abs(smooth(diff(sqrt(X_body.^2+Y_body.^2)),0.3));
                intervalMin = max(d_body)-min(d_body); %speed span
                %[xr,xy] = modeMap(d_body,5e-6);
                %[~,id]=max(xy);
                selectionTime = logical([d_body>(mode(d_body)+intervalMin*0.1); 0]);
                %plot(time_body(selectionTime),X_body(selectionTime))
               
                %Prolongate the leg to the tip
                nX_leg = X_leg;
                nY_leg = Y_leg;
                nZ_leg = Z_leg;
                for i=1:numel(X_leg)
                    rotation = transl(X_leg(i),Y_leg(i),Z_leg(i));
                    rotation(1:3,1:3) = rotz(Zr_leg(i))*roty(Yr_leg(i))*rotx(Xr_leg(i));
                    new_leg = rotation*[0 0 -Dtip 1]';
                    
                    nX_leg(i) =  new_leg(1);
                    nY_leg(i) =  new_leg(2);
                    nZ_leg(i) =  new_leg(3);
                end
                
            %Extract data
                X_data = X_leg-X_body;
                Y_data = Y_leg-Y_body;
                Z_data = Z_leg-Z_body;
                %Step frequency
%                 [X_freq, X_fft] = FFT(time_body,X_data);
%                 [Y_freq, Y_fft] = FFT(time_body,Y_data);
%                 [Z_freq, Z_fft] = FFT(time_body,Z_data);
%                 
%                 [~,id] = max(X_fft(1:end));
%                 fprintf("\tMax X %.4f ",X_freq(id+1))
%                 [~,id] = max(Y_fft(2:end));
%                 fprintf("\tMax Y %.4f ",Y_freq(id+1))
%                 [~,id] = max(Z_fft(2:end));
%                 fprintf("\tMax Z %.4f \n",Z_freq(id+1))
                
                
                %Step length
            
            %Plot data
%                 figure(1);
%                 clf(1)
%                 
%                 subplot(311); hold on;
%                 plot(time_body(selectionTime),X_data(selectionTime)*1e3)
%                 plot(time_body(selectionTime),Y_data(selectionTime)*1e3)
%                 plot(time_body(selectionTime),Z_data(selectionTime)*1e3)
%                 %plot3(X_data*1e3,Y_data*1e3,Z_data*1e3)
%                 %plot3((nX_leg-X_body)*1e3,(nY_leg-Y_body)*1e3,(nZ_leg-Z_body)*1e3)
%                 %view([1 1])
%                 %xlabel("x");ylabel("y");zlabel("z")
%                 grid;
%                 hold off
%                 subplot(312);hold on;
%                 plot(time_body,X_data*1e3)
%                 plot(time_body,Y_data*1e3)
%                 plot(time_body,Z_data*1e3)
%                 hold off;
%                 grid;
%                 subplot(313); hold on;
%                 plot(time_body(selectionTime),X_body(selectionTime)*1e3)
%                 %plot(time_body,X_body*1e3)
%                 hold off;
                %plot(power)    
                %grid;
                
                %userInput = input('Step distance :\n'); %Freq Stepdist currentStart currentStop
                         
            
           
            %Save data
                %Estimate the power of the robot
                Zpow(ld,freq,step) = mean(power(Zc1b(ld,freq,step):Zc2b(ld,freq,step)));
                
                %Estimate the speed of the robot
                %power storage
                Zspd(ld,freq,step) = sqrt((max(X_body(selectionTime))-min(X_body(selectionTime)))^2+(max(Y_body(selectionTime))-min(Y_body(selectionTime)))^2)/(max(time_body(selectionTime))-min(time_body(selectionTime)));
                %Zf(ld,freq,step) = userInput(1); %frequencies storage
                %Zd(ld,freq,step) = userInput(2); %distance of step (need the only leg graph)
                %Zs1(ld,freq,step) = userInput(1); %step start
                %Zc1(ld,freq,step) = userInput(3); %current start (compare segmented mean and global mean)
                %Zc2(ld,freq,step) = userInput(4); %current stop
        end
    end
end

%% Generate color grade 
colorGrade = zeros(numel(stepList),3);
baseColor = [21 209 87]./255;
for i=1:length(colorGrade)
   colorGrade(i,:) = baseColor.*[i/length(colorGrade) i/length(colorGrade) i/length(colorGrade)]
end

legendList = cell(1,numel(stepList));
for i=1:numel(stepList)
   legendList{i} = ['step ' num2str(stepList(i)*2) 'mm']; 
end

%% Plot power(load, freq, step)
figure; hold on;
for i=1:numel(stepList)
    surf(X(:,2:end,i),Y(:,2:end,i),Zpow(:,2:end,i),'FaceAlpha',0.5,'FaceColor',colorGrade(i,:))
    
end
legend(legendList)
grid;
hold off;

xlabel("Freq [Hz]")
ylabel("Load [g]")
zlabel("Avg. Power [W]")

%% Fit step length

ZstepFitp1 = zeros(size(X,1),size(X,2)-1);
ZstepFitp2 = ZstepFitp1;
ZstepFitp3 = ZstepFitp1;

figure(1);
for x=1:size(X,1)
   for y=2:size(X,2)
       Yfit = zeros(1,length(stepList));
       for st=1:length(stepList)
            Yfit(st) = Zpow(x,y,st)/Zpow(x,y,3);
       end
       [ft,gof] = fit(stepList'*2,Yfit','poly2');
       %fprintf('Fit %.2f \n',ft.p1)
       ZstepFitp1(x,y-1) = ft.p1;
       ZstepFitp2(x,y-1) = ft.p2;
       ZstepFitp3(x,y-1) = ft.p3;
       
%        clf;
%        hold on;
%        plot(stepList,Yfit,'or')
%        plot(ft)
%        
%        hold off;
%        pause(0.1);
%        pause(1)
   end
   
end

%Fit surfaces for prediction of A and B and C
xstData = reshape(X(:,2:end,1),[size(X,1)*(size(X,2)-1),1]);
ystData = reshape(Y(:,2:end,1)*1e-3,[size(Y,1)*(size(Y,2)-1),1]);
zstData1 = reshape(ZstepFitp1,[size(ZstepFitp1,1)*size(ZstepFitp1,2),1]);
zstData2 = reshape(ZstepFitp2,[size(ZstepFitp2,1)*size(ZstepFitp1,2),1]);
zstData3 = reshape(ZstepFitp3,[size(ZstepFitp3,1)*size(ZstepFitp1,2),1]);
ft = fittype( 'loess' );
opts = fitoptions( 'Method', 'LowessFit' );
opts.Normalize = 'on';
opts.Span = 0.95;
fitFct_steps1 = fit( [xstData, ystData], zstData1, ft ,opts);
fitFct_steps2 = fit( [xstData, ystData], zstData2, ft ,opts);
fitFct_steps3 = fit( [xstData, ystData], zstData3, ft ,opts);
plot(fitFct_steps3, [xstData, ystData], zstData3);

%Warning : fits in kg and frequency Hz

clear xstData ystData zstData1 zstData2 zstData3

%% Specific resistance (load, freq, step)
errF = 0;
errDx = 20;
figure; hold on;
for i=3:3%numel(stepList)
    V = 4*stepList(i)*(1e-3).*X_prev(:,:,i); %x4 since robot control code /2 steps
    %Zspd(:,:,i)
    %wait;
    %Vplus = 2*(2*stepList(i)+errDx)*(1e-3).*(X(:,:,i)+errF);
    %Vminus = 2*(2*stepList(i)-errDx)*(1e-3).*(X(:,:,i)-errF);
    CoT_fr = round(Zpow(:,:,i)./ ((robot_mass_fr+Y_prev(:,:,i))*1e-3*g.*V),1);  %W/mgv
    %CoTplus = Zpow(:,:,i)./ ((robot_mass+Y(:,:,i))/3*1e-3*g.*Vplus);  %W/mgv
    %CoTminus = Zpow(:,:,i)./ ((robot_mass+Y(:,:,i))/3*1e-3*g.*Vminus);  %W/mgv
    %surf(X(:,:,i),Y(:,:,i),CoTplus,'FaceAlpha',0.4,'FaceColor',[1 0 0])
    surf(V(:,:),Y_prev(:,:,i),CoT_fr,'FaceAlpha',0.5,'FaceColor',colorGrade(i,:))
    %surf(X(:,:,i),Y(:,:,i),CoTminus,'FaceAlpha',0.4,'FaceColor',[0.7 0 0])
    legend(legendList{i})
end
grid;
hold off;

xlabel("Freq [Hz]")
ylabel("Load [g]")
zlabel("CoT")


