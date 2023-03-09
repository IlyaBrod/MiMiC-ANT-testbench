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

%%
originalPower = est_POW_L6S140D05; %SELECT DESIRED POWER
originalStep = 140; %select step length
leg_nbr = 6; %select number of legs


%% Import basal power

X = zeros(1,sum(data.speedList==0)); %load data
Y = X; %power data
idx = 1;
for ispd = 1:length(data.speedList)

    if(data.speedList(ispd)==0)
        X(idx) = data.loadList(ispd);
        lcurr = data.expList(ispd).current;
        lvolt = data.expList(ispd).voltage;
        Y(idx) = mean(lcurr(:,2))*mean(lvolt(:,2));
        idx= idx+1;
    end
    

end
clear lcurr lvolt

basalFit = fit(X',Y','poly1');

fprintf("Done \n");

%Generate basal fit matrix
est_POW_basal = NaN(size(simu_Y));
for i=1:size(simu_Y,1)
    for j=1:size(simu_Y,2)
        est_POW_basal(i,j) = basalFit(simu_Y(i,j));
    end
end

surf(simu_X,simu_Y,est_POW_basal*6)%Basal power

est_POW_err = est_POW_basal*leg_nbr;

%% Step length variation error
%Estimate step error
fx = [0.4 0.8 1.2 1.5 1.9 2.3 2.7];
fy = [40 42 48 70 82 80 82];

stepLerrorFit = fit(fx',fy','smoothingspline','SmoothingParam',0.9941505464621068);
stepLerror = zeros(size(simu_X));

for i=1:size(simu_X,1)
    for j=1:size(simu_X,2)
        stepLerror(i,j) = stepLerrorFit(simu_X(i,j));
    end
end
    
%Estimate power variation due to step error
%Estimate k values for load frequency
stepPerror = zeros(size(simu_X));



for i=1:size(simu_X,1)
    for j=1:size(simu_X,2)
        fitStep1 = fitFct_steps1(simu_X(i,j),simu_Y(i,j));
        fitStep2 = fitFct_steps2(simu_X(i,j),simu_Y(i,j));
        fitStep3 = fitFct_steps3(simu_X(i,j),simu_Y(i,j));
        targetStep = originalStep - stepLerror(i,j);
        stepPerror(i,j) = originalPower(i,j)*(fitStep1*targetStep^2 + fitStep2*targetStep + fitStep3);
    end
end
est_POW_err = est_POW_err + (originalPower-stepPerror);

%% Temperature data
est_POW_err = est_POW_err + est_POW_tempVar*leg_nbr;


%% Plot comparison
figure;
hold on;
surf(simu_X,simu_Y,est_POW_err,'EdgeColor','none','FaceColor','interp')
hold on;
xlabel('X - Speed [m/s]')
ylabel('Y - Load [kg]')
zlabel('Z - Power [W]') %Cost of transport
cbar = colorbar;
ylabel(cbar, 'Estimated robot power','FontSize',12) %'Estimated robot CoT'
grid on
view(0,90)


