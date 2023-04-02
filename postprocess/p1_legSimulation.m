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

%% Define robot leg links 
clear Leg_sw

l1 = 146*1e-3;
l2 = 83*1e-3;
l3 = 53*1e-3;
l4 = 0;

tR = pi/180;
I_ax = [16149.49, 0.8, -2.55;...
          0.8, 8897.35, 136.95;...
          -2.55, 136.95, 13599.25]*1e-9; %AX18 motor inertia
I_F2 = [2946.73, 0.16, 0;...
        0.16, 3346.35, 0;...
        0, 0, 1292.01]*1e-9; %AX18 fixation inertia
mAX = 55.9e-3; %servomotor mass

L3_r = -[74.74 0 -4.06]*1e-3;
L3_I = [16835.71, 0.14, 0;...
          0.14, 95062.05, 0;...
          0, 0, 84081.04] * (1e-9);
L2_r =  -(([12.72 0 1.02]*55.9 + [(36*2+9.2-12.72) 0 1.02]*55.9)/(55.9*2))*1e-3;
L2_I = ihuygens(I_ax,([12.72 0 1.02]*1e-3-L2_r),mAX)+ihuygens(I_ax,([(36*2+9.2-12.72) 0 1.02]*1e-3-L2_r),mAX);      

L1_r = -((7.49*[15.40 0 0]+7.49*[(26.5*2-15.40) 0 0])/(2*7.49))*1e-3;
L1_I = ihuygens(I_F2,([15.40 0 0]*1e-3-L1_r),7.49e-3)+ihuygens(I_F2,([(26.5*2-15.40) 0 0]*1e-3-L1_r),7.49e-3);

L1_m = 7.49*2*1e-3;
L2_m = mAX*2;
L3_m = 39.23*1e-3;

% Links definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SWING
L_sw(1) = Revolute('d', 0, ...   % link length (Dennavit-Hartenberg notation)
'a', 53*1e-3, ...               % link offset (Dennavit-Hartenberg notation)
'alpha', pi/2, ...        % link twist (Dennavit-Hartenberg notation)
'I', [L1_I(1,1), L1_I(2,2), L1_I(3,3), L1_I(1,2), L1_I(2,3), L1_I(1,3)], ... % inertia tensor of link with respect to center of mass I = [L_xx, L_yy, L_zz, L_xy, L_yz, L_xz]
'r', L1_r, ...       % distance of ith origin to center of mass [x,y,z] in link reference frame
'm', L1_m);   
   
L_sw(2) = Revolute('d', 0, ...   % link length (Dennavit-Hartenberg notation)
'a', 83*1e-3, ...               % link offset (Dennavit-Hartenberg notation)
'alpha', 0, ...        % link twist (Dennavit-Hartenberg notation)
'I', [L2_I(1,1), L2_I(2,2), L2_I(3,3), L2_I(1,2), L2_I(2,3), L2_I(1,3)], ... % inertia tensor of link with respect to center of mass I = [L_xx, L_yy, L_zz, L_xy, L_yz, L_xz]
'r', L2_r, ...       % distance of ith origin to center of mass [x,y,z] in link reference frame
'm', L2_m);

L_sw(3) = Revolute('d', 0, ...   % link length (Dennavit-Hartenberg notation)
'a', 146*1e-3, ...               % link offset (Dennavit-Hartenberg notation)
'alpha', 0, ...        % link twist (Dennavit-Hartenberg notation)
'I', [L3_I(1,1), L3_I(2,2), L3_I(3,3), L3_I(1,2), L3_I(2,3), L3_I(1,3)], ... % inertia tensor of link with respect to center of mass I = [L_xx, L_yy, L_zz, L_xy, L_yz, L_xz]
'r', L3_r, ...       % distance of ith origin to center of mass [x,y,z] in link reference frame
'm', L3_m); 

figure(1);
Leg_sw=SerialLink(L_sw,'name','AXLeg');
%Leg_st.plot([0 45 -45 0]*tR,'arrow','wrist','base','jaxes')
Leg_sw.plot([0 0 -90]*tR,'arrow','wrist','base','jvec')
clear L_sw

%% Power estimation

RES = 50; %Power param resolution
simu_speeds = 0.1:(0.5/RES):0.6;
simu_freqs = 0.5:2.5/RES:3;
simu_loads = 0:(0.5/RES):0.5;
simu_robot_mass = 2.183;
duty_factor = 0.5;
GaitList;gaitSelect = gait_wave;COT_INVERT_GAIT = 0; 
legs_nbr = numel(gaitSelect);

stepH = -40e-3; %negative offset compared to body high
stepL = 100e-3; %centered on the middle axe of the leg [100 default]
stepD = 180e-3; %distance from body to the leg tip
simu_res = 100; %Time resolution

ptz_t = [-0.2 0  0.2 0.3 0.4 0.45 0.6  0.8  1];
ptz_v = [0 1  1    1   1   1   0.5   0    1];

pty_t = [-0.15 0 0.25 0.5 0.75 1.0];
pty_v = [0 -0.5 0 0.5 0 -0.5];

interp_time = 0:1/(simu_res):1;
interpZ_value = spline(ptz_t,ptz_v,interp_time);
interpY_value = spline(pty_t,pty_v,interp_time);
interpX_value = stepD*ones(size(interpZ_value));

X_simu = interpX_value';
Y_simu = interpY_value'*stepL;
Z_simu = interpZ_value'*stepH;   

X_simu = [X_simu ; X_simu(2:end) ; X_simu(2:end)];
Y_simu = movingAverage([Y_simu ; Y_simu(2:end) ; Y_simu(2:end)],10);
Z_simu = movingAverage([Z_simu ; Z_simu(2:end) ; Z_simu(2:end)],10);

segment = round(length(Z_simu)/3);
Z_simu = Z_simu(segment:segment*2);
Y_simu = Y_simu(segment:segment*2);
X_simu = X_simu(segment:segment*2);
figure;hold on;plot(interp_time,Z_simu);plot(interp_time,Y_simu);hold off;
figure;plot(Y_simu,Z_simu);

tselect_sw = (Z_simu > min(Z_simu)*0.95); %swing
tselect_st = ~tselect_sw; %stance
%tselect_sw = dilateMask(tselect_sw,4,'Cycle');
%tselect_st = dilateMask(tselect_st,4,'Cycle');

figure(1);
clf;
hold on;
plot3(X_simu(tselect_sw)*1e3,Y_simu(tselect_sw)*1e3,Z_simu(tselect_sw)*1e3,'o','LineWidth',2)
plot3(X_simu(tselect_st)*1e3,Y_simu(tselect_st)*1e3,Z_simu(tselect_st)*1e3,'x','LineWidth',2)
hold off
legend('swing phase','stance phase')
grid;
view(90,0)

[X,Y] = meshgrid(simu_freqs,simu_loads); %simu_speeds
simu_POW = zeros(size(X));
simu_POWe = simu_POW;
%IK
q_swO = Leg_sw.ikine(transl([X_simu Y_simu Z_simu-75e-3]),[],'mask',[1 1 1 0 0 0]);
%q_st = Leg_st.ikine(transl([X_simu Y_simu -(Z_simu-75e-3)]),[],'mask',[1 1 1 0 0 1]);
%figure; hold on;plot(scaledTime(tselect),q_swO(tselect,:),'ob'); plot(scaledTime(tselect_inv),q_swO(tselect_inv,:),'or','LineWidth',2);hold off;
[T_sw] = Leg_sw.fkine(q_swO);


tic
pBar = progressBar([numel(simu_loads) numel(simu_speeds)],5);
for spdi=1:numel(simu_freqs)
    for ldi=1:numel(simu_loads)
        pBar.update
        sp = simu_speeds(spdi);
        ld = simu_loads(ldi);
        
        %define trajectory & build time vector
            freq = simu_freqs(spdi); %sp/(2*stepL); %walking cycle frequency [Hz]
            Ts = 1/freq;
            Tech_st = duty_factor*Ts/(sum(tselect_st));
            Tech_sw = (1-duty_factor)*Ts/(sum(tselect_sw)-1);
            scaledTime = zeros(1,numel(interp_time)+4);        

            for t=1:numel(scaledTime)
               if(tselect_st(mod(t-2-1,numel(tselect_st))+1))
                  scaledTime(t) = scaledTime(mod(t-1-1,numel(scaledTime))+1) + Tech_st; 
               else
                   scaledTime(t) = scaledTime(mod(t-1-1,numel(scaledTime))+1) + Tech_sw;
               end
            end
            
            %scaledTime(end) %should equals Ts scaledTime(end-2)-scaledTime(3)  
        
        %Calculate velocity & acceleration
            q = [q_swO(end-2:end-1,:) ; q_swO ; q_swO(2:3,:)]; %extend q_swO by loop
            [qd,q_time] = center_derivate(q,scaledTime');
            [qdd,q_time] = center_derivate(medfilt1(qd),q_time);
            qdd = medfilt1(qdd);
        
        %Generate load pattern
        Li = zeros(legs_nbr,simu_res+1);
        idx = find(tselect_st); %step offset due to path generation t(idx(1))

        for i=1:legs_nbr
            Theta = gaitSelect(i);
            Li(i,:) = 1-((vz_function_rectWave((q_time-Theta*Ts-q_time(idx(1)))*pi/Ts) ...
                    .*vz_function_rectWave((q_time-Theta*Ts-duty_factor*Ts-q_time(idx(1)))*pi/Ts)+1)/2);
        end
        L = medfilt1(sum(Li));
%         figure ; hold on;
%         plot(L); plot(tselect_st);
%         hold off;

        %Calculating each different type of power value
        loadDiv_st = unique(L);
        loadDiv_st_count = numel(loadDiv_st);
            
        %Dynamics
        %Swing phase dynamics
        DynQ_sw = zeros(size(qdd));

        parfor i=1:length(qdd)
            crrP = q(2+i,:);
            crrDP = qd(1+i,:);
            crrDDP = qdd(i,:);
            DynQ_sw(i,:) = DynQ_sw(i,:)+Leg_sw.gravload(crrP);        
            DynQ_sw(i,:) = DynQ_sw(i,:)+transpose(Leg_sw.inertia(crrP)*crrDDP');
            DynQ_sw(i,:) = DynQ_sw(i,:)+ transpose(Leg_sw.coriolis(crrP,crrDP)*crrDP');
        end

        %Stance phase dynamics
        % Supported phase = aerial + forces actions
        DynQ_st = DynQ_sw;

        for i=1:length(qdd)
               crrP = q(2+i,:);
               crrDP = qd(1+i,:);
               crrDDP = qdd(i,:);
               
               Li = L(i);
               if(Li==0), continue
               else, payload = ((legs_nbr-Li)*data.mass+(ld+0.38))/Li; %(3*data.mass+(ld+0.38))/3;
               end
               FextLD = GRAVITY*payload;
               
               %calculation of center of mass
               tM = vz_robot_getTransformMat(Leg_sw,q_swO(i,:),'end',1);
               Pt = tM * [zeros(4,3) [-L1_r 1]'];
               Pt1 = Pt(1:3,end);
               %plot3(Pt1(1),Pt1(2),Pt1(3),'or','MarkerFaceColor','blue');
               tM = vz_robot_getTransformMat(Leg_sw,q_swO(i,:),'end',2);
               Pt = tM * [zeros(4,3) [-L2_r 1]'];
               Pt2 = Pt(1:3,end);
               %plot3(Pt2(1),Pt2(2),Pt2(3),'or','MarkerFaceColor','blue');
               tM = vz_robot_getTransformMat(Leg_sw,q_swO(i,:),'end',3);
               Pt = tM * [zeros(4,3) [-L3_r 1]'];
               Pt3 = Pt(1:3,end);
               %plot3(Pt3(1),Pt3(2),Pt3(3),'or','MarkerFaceColor','blue');
               XYZ = tM(1:3,end);%leg tip position
               COM = vz_phys_centerOfMass([Pt1' ; Pt2' ; Pt3'; [0 0 0]],[L1_m;L2_m;L3_m;payload]);

               %XYZ2 = T_sw(i).tv;
               %alpha = atan2(abs(XYZ2(3)),abs(XYZ2(1)));
               %Ft2 = -FextLD*cos(alpha)/sin(alpha);
               %Change to the leg tip frame and compute force
               Ft = FextLD*(XYZ(1)-COM(1))/(XYZ(3));
               %disp(Ft2-Ft)
               %Rotating to end-effector frame
               tMat = vz_robot_getTransformMat(Leg_sw,q_swO(i,:)); %getting rot matrix
               tMat = tMat(1:3,1:3);
               Fcurr = tMat\[Ft 0 FextLD]'; %we need to inverse the rot matrix
               %DynQ_st(i,:) = DynQ_st(i,:)+ Leg_sw.rne(q(i+2,:),qd(i+1,:),qdd(i,:),'fext',[Fcurr(1:3,1)' 0 0 0]);
               DynQ_st(i,:) = DynQ_st(i,:)+ transpose(Leg_sw.jacobe(crrP)'*[Fcurr(1:3,1)' 0 0 0]');
        end
        
        %torque to current (gear ratio included)
        %Nm2A = 2.2/1.83*12;%Newton-meter to Amps
        %SimuCurr_st = sum(abs(DynQ_st),2)*Nm2A.*tselect_st;
        %SimuCurr_sw = sum(abs(DynQ_sw),2)*Nm2A.*tselect_sw;
        %SimuCurrTotal_elec = SimuCurr_st+SimuCurr_sw;
        
        %Total robot power calculation
        SimuCurr_st = sum(abs(DynQ_st.*qd(2:end-1,:)),2).*tselect_st;
        SimuCurr_sw = sum(abs(DynQ_sw.*qd(2:end-1,:)),2).*tselect_sw;
        SimuCurrTotal = SimuCurr_st+SimuCurr_sw; % + basalFit(payload)
        
        RobotPower = Tech_sw*sum(SimuCurrTotal(tselect_sw)) + ...
               Tech_st * sum(SimuCurrTotal(tselect_st));
        RobotPower = legs_nbr*RobotPower/Ts;
        
%       SimuCurrTotal = zeros(size(tselect(3:end-2)));
%       SimuCurrTotal(tselect(3:end-2)) = SimuCurr_sw(tselect(3:end-2)');
%       SimuCurrTotal(tselect_inv(3:end-2)) = SimuCurr_st(tselect_inv(3:end-2));
        %SimuCurr_sw = abs(DynQ_sw*Nm2A);
        %SimuCurrTotal_sw = SimuCurr_sw(:,1) + SimuCurr_sw(:,2) + SimuCurr_sw(:,3);

        %SimuCurr_st = abs(DynQ_st*Nm2A);
        %SimuCurrTotal_st = SimuCurr_st(:,2) + SimuCurr_st(:,3);% + SimuCurr_st(:,3); SimuCurr_st(:,1)

        %SimuCurrTotal = zeros(size(scaledTime));
        %SimuCurrTotal(tselect(3:end-2)) = SimuCurrTotal_sw';
        %SimuCurrTotal(tselect_inv(3:end-2)) = SimuCurrTotal_st';
        
        %figure;hold on;plot(SimuCurrTotal);plot(SimuCurr_st,'--');plot(SimuCurr_sw,'--');plot(tselect_inv(3:end-2)) ;hold off;
        
        %RobotPower = legs_nbr*mean(SimuCurrTotal);
        simu_POW(ldi,spdi) = RobotPower;
        %simu_POWe(ldi,spdi) = legs_nbr*mean(SimuCurrTotal_elec);
        
    end
end
disp(toc/60)
%figure; plot([SimuCurrTotal(tselect(3:end-2)) SimuCurrTotal(tselect_inv(3:end-2))])
simu_X = X; simu_Y = Y;

figure;
surf(simu_X,simu_Y,simu_POW,'EdgeColor','none','FaceColor','interp')
xlabel('X - Speed [m/s]')
ylabel('Y - Load [kg]')
zlabel('Z - Robot power [W]')
cbar = colorbar;
ylabel(cbar, 'Simulated robot power','FontSize',12)
grid on
view(0,90)
title(["Simulated robot power"])
