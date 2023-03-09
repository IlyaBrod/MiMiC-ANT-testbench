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
SETTINGS

%% CONVERTING UNITS to [kg - sec - N]
expNbr = length(data.expList);
data.loadList = round(data.loadList,MASS_DECIMALS) * 1e-3;
data.mass = data.mass * 1e-3;
data.timeList = data.timeList * 60;


%% LEG TIP ESTIMATION
% Prolongate the tip of the leg caused by the impossibility to put
% a tracking marker on the tip of the leg

barFish = progressBar(expNbr,10);
for idExp=1:expNbr
    barFish.update;
    
    currExp = data.expList(idExp);
    %Leg tip estimation
    lMotion = currExp.legMotion;
    mdataSize = size(lMotion,1);
    tipMotion = zeros(mdataSize,3);
    for idPos=1:mdataSize
        pt1 = lMotion(idPos,2:4);
        pt2 = lMotion(idPos,5:7);
        pt3 = lMotion(idPos,8:10);
        tip = markerTriangleOffset(pt1,pt2,pt3,LEG_TIBIA_OFFSET);
        tipMotion(idPos,:) = tip;
    end
    currExp.legMotion_ = [lMotion(:,1) tipMotion lMotion(:,2:end)];
end


%% LEG FORCE ESTIMATION
fprintf("Load calibration file\n")
fPlate = forcePlate;
fPlate.init_calibrationZ(calib)
barFish = progressBar(expNbr,5);
for idExp=1:expNbr
    barFish.update;
    currExp = data.expList(idExp);
    fPlate.Dx = currExp.calib_dx;
    fPlate.Dy = currExp.calib_dy;
    fPlate.Dz = currExp.calib_dz;

    
    force = downsample(currExp.force,100);
    forceNorm = zeros(size(force));
    forceNorm(:,1) = force(:,1);
    lMotion = currExp.legMotion_;
    
    %Set same sampling rate
    [~, lMotion] = resampleData(lMotion, force);
    
    for t=1:length(force)
        x = lMotion(t,1);
        y = lMotion(t,2);
        %size(lMotion)
        voltage = force(t,2);
        forceNorm(t,2) = fPlate.getForce_norm(x,y,voltage);
    end
    currExp.force_ = forceNorm;
end

%% APPLY MODIFICATIONS
for i=1:expNbr
   data.expList(i).format = {'SI'}; 
end

dataBase.applyChanges(data);
