%% MIT License
% 
% Copyright (C) 2022  Ilya Brodoline
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% 
% contact: ilya.brodoline@univ-amu.fr
%%

classdef MiMicAntBench < handle
    %MiMicAntBench Manage the test bench sensors and output
    %                      Data storage format:
    %scenario X data: 
    %scxData =             name
    %                      mass
    %                      tempList = [exp1temp, exp2temp...]
    %                      loadList = [exp1load, exp1load...]
    %                      speedList = [speed1,speed2]
    %                      expList = [exp1Container, exp2Container...]
    %
    %
    %Experiment Y of scenario X data:
    %expxContainer =             expName
    %                            id (same for a speed range)
    %                            speed
    %                            force = [Yforce vector , Zforce vector]
    %                            temp = [Temp vector]
    %                            voltage = [Voltage vector]
    %                            current = [Current vector]
    %                            legMotion = [Leg tip XYZ , ... , Leg base XYZ]
    %                            groundMotion = [XYZ point Y increase & X increase order]
    %                            date = "15-Mar-2022 14:35:38"
    
    properties
        %Test bench modules
        FLIR;
        ADC;
        controller;
        resultTab;
        mCAP;
        fPlate;
        
        %Experimental data
        expAdd=0;
        sc1Data;
        sc2Data;
        sc3Data;

        %Debug calibration
        debugData = containers.Map({'layers','weight_matrix','layers_X','layers_Y'},{0,0,0,0},'UniformValues',false)
    end
    
    properties (Access = private)
        alert;
    end
    
    methods
        function obj = MiMicAntBench()
            %MiMicAntBench Init the test bench sensors            
            SETTINGS
            
            obj.sc1Data = scData("Walking free");            
            obj.sc2Data = scData("Walking constrain");
            obj.sc3Data = scData("Lifting");
            
            obj.sc1Data.mass = LEG_MASS;
            obj.sc2Data.mass = LEG_MASS;
            obj.sc3Data.mass = LEG_MASS;
            
            fprintf("INIT sensors...\n");
            fprintf("\tFLIR thermal module\n")
            obj.FLIR = thermalCamera(thermal_PORT);

            fprintf("\tZ Force & Current sensors\n")
            obj.ADC = NationalInstrument(ni_RECORD_RATE,ni_RECORD_TIME);
            obj.ADC.addChannel('Current',ni_current_PORT);
            obj.ADC.addChannel('Zforce',ni_sensZ_PORT);
            %obj.ADC.addChannel('Voltage',ni_voltage_PORT);
            obj.fPlate = forcePlate();
            obj.fPlate.init_calibrationZ(fPlate_calib);
            
            if(fPlate_STIFFNESS>0), obj.fPlate.stiffness=fPlate_STIFFNESS;end
            
            fprintf("\tLeg controller\n")
            obj.controller = robotController("COM5",9600,0); %seriallist
            
            fprintf("\tMotion capture\n")
            obj.mCAP = MotionCapture;
            
            obj.resultTab = resultTable();
            obj.alert = alertBox;
            obj.alert.sound("bip");
            fprintf("Initialisation of the test bench DONE\n");
        end

        function saveRAW(obj,name)
            SETTINGS
            
            sc1Data = obj.sc1Data;
            save([CACHE_EXPORT '/' name '_sc1.mat'],'sc1Data')
            delete(sc1Data);
            fprintf("25%%\n");
            sc2Data = obj.sc2Data;
            save([CACHE_EXPORT '/' name '_sc2.mat'],'sc2Data')
            delete(sc2Data);
            fprintf("50%%\n")
            sc3Data = obj.sc3Data;
            save([CACHE_EXPORT '/' name '_sc3.mat'],'sc3Data')
            delete(sc3Data);
            fprintf("75%%\n")
            debugData = obj.debugData;
            save([CACHE_EXPORT '/' name '_debug.mat'],'debugData')
            delete(debugData);
            fprintf("100%%\n")
            fprintf("RAW data saved to cache folder\n");
        end
        
        function addExperiment(obj)
            obj.expAdd = 1;
        end
        
        function cancelAddExperiment(obj)
           obj.expAdd = 0; 
        end
        
        function resetResults(obj,confirm)
            %resetResults Construct and initialize the result table
            % Write confirm=1 to proceed.
            if(confirm == 1)
                obj.resultTab = resultTable();
            end
        end
        
        function force = measureZForce(obj,markerID)
            data = obj.ADC.read(); %record adc
            idx = obj.ADC.getIndex('Zforce');
            voltage = data(idx);
            [x,y,~] = obj.mCAP.getPoint(markerID); %get x,y
            force = obj.fPlate.getForce_norm(x,y,voltage); %calculate force
        end
        
        function run(obj,scenario)
            %run Run a specific scenario to estimate specific results.
            %   argument scenario number:
            %(-2)- Force Plate Calibration Horizontal
            %(-1)- Force Plate Calibration Vertical
            % 1- Walk 
            % (Stride length ; Stride max freq ; Innefective stance ;
            % step accuracy ; CoT ; Temp inc.)
            % 2- Walk with resistance
            % (Max. power ; Max tan ; Temp inc. ;
            % 3- Weight lifting 
            % (Max norm stress ; Temp. inc.)

            SETTINGS
            switch scenario
                    
                case -2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CALIBRATION Y
                    
                    fprintf(obj.alert.msg_calibration2_warn);
                    wait; obj.alert.sound("message");
                    [~,y11,~] = obj.mCAP.getPoint(mCAP_ID('TD_1')); %Getting threadmill 
                    [~,y21,~] = obj.mCAP.getPoint(mCAP_ID('TD_2')); %markers position
                    [~,y31,~] = obj.mCAP.getPoint(mCAP_ID('TD_3'));
                    fprintf("Pull back the dynamometer (at least 5cm)\n");
                    F = input('Enter measured value [N]: ');
                    [~,y12,~] = obj.mCAP.getPoint(mCAP_ID('TD_1')); %Getting threadmill 
                    [~,y22,~] = obj.mCAP.getPoint(mCAP_ID('TD_2')); %markers position
                    [~,y32,~] = obj.mCAP.getPoint(mCAP_ID('TD_3'));
                    displace = mean(abs([y11-y12 y21-y22 y31-y32])); %Mean displace
                    K = single(F)/(displace*1e-3);
                    obj.fPlate.stiffness = K;
                    %do same
                    fprintf("Threadmill stiffness coefficient k=%3.4f\n",K)
                    fprintf("Calibration DONE\n");wait;
                    
                    
                case -1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CALIBRATION Z
                    
                    fprintf(obj.alert.msg_calibration_warn);
                    wait; obj.alert.sound("message");
                    fprintf(obj.alert.msg_calibration_init);
                    obj.alert.sound("message"); wait;
                    fprintf("\n\nSTARTING\n");
                    obj.alert.sound("start")
                    
                    %Start
                    %Dimension detection
                    fprintf(obj.alert.msg_calibration_corner1);
                    obj.alert.sound("message"); wait;
                    obj.mCAP.updateLabels;
                    [x1,y1,~] = obj.mCAP.getPoint(mCAP_ID('MK_CALIB_2'));
                    [x2,y2,~] = obj.mCAP.getPoint(mCAP_ID('MK_CALIB_3'));
                    [x3,y3,~] = obj.mCAP.getPoint(mCAP_ID('MK_CALIB_4'));
                    fPlate_calib.xOffset = min([x1 x2 x3]); % TODO add NaN detection
                    fPlate_calib.yOffset = min([y1 y2 y3]);
                    fPlate_calib.xSize = max([x1 x2 x3]) - fPlate_calib.xOffset;
                    fPlate_calib.ySize = max([y1 y2 y3]) - fPlate_calib.yOffset;
                    obj.fPlate.init_calibrationZ(fPlate_calib);
                    fprintf("Calculated plate dimensions: X:%3.2fmm Y:%3.2fmm\n",obj.fPlate.calibration.xSize,obj.fPlate.calibration.ySize);
                    obj.alert.sound("message")
                    fprintf(obj.alert.msg_calibration_corner2); wait;
                    
                    %0.0gramms calibration measure
                    layer_count = 1;
                    weight_matrix = zeros(fPlate_calib_minimal_layers+1,1);
                    layers = zeros(fPlate_calib.Xresolution,fPlate_calib.Yresolution,fPlate_calib_minimal_layers+1);
                    layers_X =zeros(fPlate_calib.Xresolution,fPlate_calib.Yresolution,fPlate_calib_minimal_layers+1);
                    layers_Y =zeros(fPlate_calib.Xresolution,fPlate_calib.Yresolution,fPlate_calib_minimal_layers+1);
                    
                    fprintf(obj.alert.msg_calibration_zero); wait;
                    fprintf("Calibrating ...");
                    idx = obj.ADC.getIndex('Zforce')+1;
                    data = obj.ADC.recordSeconds(5);
                    layers(:,:,1) = ones(fPlate_calib.Xresolution,fPlate_calib.Yresolution);
                    layers(:,:,1) = layers(:,:,1).*mean(data(:,idx));
                    for i=1:fPlate_calib.Xresolution %save 0g locations
                        for j=1:fPlate_calib.Yresolution
                            [x,y] = obj.fPlate.reverseIndex(i,j);
                            layers_X(i,j,1) = x;
                            layers_Y(i,j,1) = y;
                        end
                    end
                    layer_count = layer_count +1;
                    fprintf("DONE\n");
                    obj.alert.sound("bip");
                    
                    % Other weights measurment
                    markerID = input('Enter weight marker ID (+ Labeled - Unlabeled): ');
                    fig = figure(1); %Preview figure
                    while(layer_count<=fPlate_calib_minimal_layers+1)
                        layer = zeros(fPlate_calib.Xresolution,fPlate_calib.Yresolution);
                        fprintf(obj.alert.msg_calibration_addmass);
                        fprintf("<Warning> Do not close the figure\n");
                        curr_weight = input('Mass value [g]: ');
                        weight_matrix(layer_count) = curr_weight*1e-3;
                        obj.alert.sound("message")
                        fprintf(obj.alert.msg_calibration_process);
                        calib_callback_state = false;
                        surf(abs(layer./(max(max(layer))+0.001)),'EdgeColor','interp','FaceAlpha',0.8,'LineWidth',20); %Display calibration state
                        view(90,-90)
                        while(calib_callback_state == false)
                            %Get DATA
                            [x,y,~] = obj.mCAP.getPoint(markerID);
                            if(isnan(x) || isnan(y)), continue, end
                            HQ = 1;
                            while(HQ>CALIB_VAR) %Automatic detection of weight placement
                                data = obj.ADC.recordSeconds(CALIB_DELAY);
                                HQ = std(data(:,idx));
                            end
                            %Update LAYER
                            [i, j] = obj.fPlate.getIndex(x, y);
                            if(i<=0 || i> size(layer)*[1;0]  || j<=0 || j>size(layer)*[0;1])
                                fprintf("<Marker out of plate %d %d>\n",i,j);
                            else
                                if(layer(i,j)>0) %If the value already exists, we do a mean
                                    layer(i,j) = mean([layer(i,j) mean(data(:,idx))]);
                                    layers_X(i,j,layer_count) = mean([layers_X(i,j,layer_count) x]);
                                    layers_Y(i,j,layer_count) = mean([layers_Y(i,j,layer_count) y]);
                                else %Set the first value
                                    layer(i,j) = mean(data(:,idx));
                                    layers_X(i,j,layer_count) = x;
                                    layers_Y(i,j,layer_count) = y;
                                end
                                %Plot PREVIEW
                                clf(fig) 
                                hold on;
                                surf(abs(layer./(max(max(layer))+0.001)),'EdgeColor','interp','FaceAlpha',0.8,'FaceColor','interp','LineWidth',20); %Display calibration state
                                plot3(j,i,max(max(layer))*2,'-o','Color','r','MarkerSize',10, 'MarkerFaceColor',[128 98 245]./255)
                                view(90,-90)
                                hold off;
                            end
                            if(sum(sum(layer==0))==0) %Check if layer is full
                                calib_callback_state = true;
                                obj.alert.sound("bip");
                            end
                        end
                        %SAVE calibration layer
                        layers(:,:,layer_count) = layer;
                        layer_count = layer_count+1;
                        surf(abs(layer./(max(max(layer))+0.001)),'EdgeColor','interp','FaceAlpha',0.8,'FaceColor','interp','LineWidth',20); %Display calibration state
                    end
                    surf(abs(layer./(max(max(layer))+0.001)),'EdgeColor','interp','FaceAlpha',0.8,'FaceColor','interp','LineWidth',20); %Display calibration state
                    view(90,-90)
                    fprintf("Measurment finished\n"); wait;
                    
                    %CALIBRATION: multivariable regression 
                    % Y = X.B B=[b1, b2, b3, b4] X = [1 volt X Y]
                    % Estimation is given by calculating B
                    
                    nsize = fPlate_calib_minimal_layers+1;
                    X = zeros(fPlate_calib.Xresolution*fPlate_calib.Yresolution*nsize,4);
                    Y = zeros(fPlate_calib.Xresolution*fPlate_calib.Yresolution*nsize,1);
                    idx =1; %Index of column
                    
                    for ly=1:nsize %Constructing X matrix
                        curr_weight = weight_matrix(ly);
                        curr_layer = layers(:,:,ly);
                        for i=1:fPlate_calib.Xresolution
                            for j=1:fPlate_calib.Yresolution
                                Y(idx,1) = curr_weight;
                                X(idx,:) = [1 curr_layer(i,j) layers_X(i,j) layers_Y(i,j)];
                                idx = idx+1;
                            end
                       end
                    end
                    %Solving
                    [B,~,r,rint,~]  = regress(Y,X);
                    %Display quality of calibration
                    contain0 = (rint(:,1)<0 & rint(:,2)>0);
                    idx = find(contain0==false); %Contains outliers indexes
                    clf(fig)
                    subplot(211)
                    hold on
                    scatter(Y,r)
                    scatter(Y(idx),r(idx),'b','filled')
                    xlabel("Mass")
                    ylabel("Residuals")
                    title("Regression residuals")
                    hold off
                    subplot(212)
                    hold on
                    scatter3(X(:,3),X(:,4),Y,'filled','MarkerFaceColor',[0 .75 .75])
                    scatter3(X(:,3).*contain0,X(:,4).*contain0,Y.*contain0,'filled','MarkerFaceColor',[0.75 0.2 0])
                    view(-30,10)
                    hold off
                    xlabel("X [mm]")
                    ylabel("Y [mm]")
                    zlabel("Weight [g]")
                    title("Residuals position")
                    %Save calibration
                    fPlate_calib.B = B;
                    fPlate_calib.resetCalibration = false;
                    obj.fPlate.init_calibrationZ(fPlate_calib);
                    %Save debug data
                    obj.debugData('layers') = layers;
                    obj.debugData('layers_X') = layers_X;
                    obj.debugData('layers_Y') = layers_Y;
                    obj.debugData('weight_matrix') = weight_matrix;
                    fprintf("Calibration DONE\n");
                    obj.alert.sound("taskend"); wait;
                    close(1)
                    
                case 1   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SCENARIO 1
                    
                    fprintf(obj.alert.msg_scenario_1_warn);
                    obj.alert.sound("message"); wait;
                    
                    %Init vars
                    %exp counter increase
                     if(obj.expAdd~=1)
                            expID = obj.sc1Data.id_count +1;
                            obj.sc1Data.id_count = obj.sc1Data.id_count +1;
                     else
                            expID = obj.sc1Data.id_count;
                     end
                     
                    forceIdx = obj.ADC.getIndex('Zforce');
                    currIdx = obj.ADC.getIndex('Current');
                    %voltIdx = obj.ADC.getIndex('Voltage');
                    obj.ADC.setRecordTime(SC1_TIME*60);
                    
                    %Data size definition (1 mCap, 2 thermal,3 ADC)
                    sSize_thread1 = round(SC1_TIME*60*1.1/THREAD1_TECH); %Number of iterations
                    sSize_thread2 = round(SC1_TIME*60/THREAD2_TECH); %Number of iterations
                    speedSequence = SC1_MIN_SPEED:(SC1_MAX_SPEED-SC1_MIN_SPEED)/(SC1_SPEED_STEPS-1):SC1_MAX_SPEED;
                    th1_dly = THREAD1_TECH*1000;
                    th2_dly = THREAD2_TECH*1000;
                    th1_max = SC1_TIME*60;
                    th2_max = SC1_TIME*60;
                    th3_max = SC1_TIME*60;
                    fprintf("Speed sequence :");
                    disp(speedSequence);
                    
                    %Start process
                    obj.alert.sound("start")
                    delete(obj.FLIR);
                    delete(obj.mCAP);
                    poolObj = gcp('nocreate'); delete(poolObj);
                    poolObj = parpool(2,'AttachedFiles',{'SETTINGS.m'});
                    
                    %obj.FLIR.lockUntil_max(ROOM_TEMP,ROOM_TEMP_ERR); %Wait motor cooldown
                    %obj.FLIR.lockUntil_mean(ROOM_TEMP,ROOM_TEMP_ERR);%Wait room temp
                    
                    for iSpd = 1:SC1_SPEED_STEPS
                        obj.alert.sound("bip")
                        
                        %Init measure data
                        data_mocap = zeros(sSize_thread1,4,2);
                        data_temp = zeros(sSize_thread2,2,2);
                        
                        %Set the current speed
                        currSpeed = speedSequence(iSpd);
                        obj.controller.setSpeed(currSpeed);

                        fprintf("Scenario %d\n",scenario);
                        fprintf("Walk speed %d/%d \n",iSpd,SC1_SPEED_STEPS);
                        fprintf("Ambiant temps %d °C\n",ROOM_TEMP);
                        fprintf("Load %d g\n",LEG_MASS+LEG_LOAD);
                        fprintf("Processing ... \n");
                        tic %start timer
                        obj.controller.start(1);
                        
                        capID = mCAP_ID('TB_1');
                        camPort = thermal_PORT;
                        
                        %Start
                        obj.ADC.recordThread(th3_max);
                        parfor thread=1:2
                           fprintf("\tLocal thread %d started\n",thread);
                           tic
                           switch(thread)                 
                               case 1
                                   %mocap
                                    localmCAP = MotionCapture;
                                    timePrev = toc;
                                    for t=1:sSize_thread1
                                        [x,y,z] = localmCAP.getPoint(capID);
                                        time = toc;
                                        data_mocap(t,:,thread) = [time x y z];
                                        delay(th1_dly-(toc-timePrev)*1000);
                                        timePrev = toc;
                                        %if(toc>th1_max), break;end
                                    end
                                   
                               case 2
                                   localCam = thermalCamera(camPort);
                                   localCam.start
                                   timePrev = toc;
                                   for tt=1:sSize_thread2
                                        time = toc;
                                        localCam.capture;
                                        data_temp(tt,:,thread) = [time max(max(localCam.frame))];
                                        delay(th2_dly-(toc-timePrev)*1000);
                                        timePrev = toc;
                                        %if(toc>th2_max), break;end
                                   end
                                   
                           end
                            fprintf("\tLocal thread %d terminated (%4.4f sec)\n",thread,toc);
                        end
                        
                        wait(obj.ADC.DAC);
                        fprintf("\tADC terminated\n");
                        obj.controller.stop;
                        fprintf("Done\n");
                        fprintf("Step time %4.4f\n",toc);
                        
                        %Save data block
                        obj.sc1Data.loadList = [obj.sc1Data.loadList (LEG_MASS+LEG_LOAD)];
                        obj.sc1Data.tempList = [obj.sc1Data.tempList ROOM_TEMP];
                        obj.sc1Data.speedList = [obj.sc1Data.speedList currSpeed];
                        obj.sc1Data.timeList = [obj.sc1Data.timeList SC1_TIME];
                        
                        expContainer = expData(strcat("SC1_L",int2str(LEG_MASS+LEG_LOAD),"_T",int2str(ROOM_TEMP),"_S",int2str(currSpeed)));
                        expContainer.date = datetime;
                        expContainer.id = expID;
                        expContainer.speed = currSpeed;
                        expContainer.force = [obj.ADC.time(:,1) ,obj.ADC.data(:,forceIdx)];
                        expContainer.temp = data_temp(:,:,2);
                        expContainer.voltage = motor_VOLTAGE;%[obj.ADC.time(:,1) ,obj.ADC.data(:,voltIdx)];
                        expContainer.current = [obj.ADC.time(:,1) ,obj.ADC.data(:,currIdx)];
                        expContainer.legMotion = data_mocap(:,:,1);
                        expContainer.groundMotion = 0;
                        

                        obj.ADC.stopThread;
                        
                        %obj.FLIR.lockUntil_max(ROOM_TEMP,ROOM_TEMP_ERR); %Wait motor cooldown
                        %obj.FLIR.lockUntil_mean(ROOM_TEMP,ROOM_TEMP_ERR);%Wait room temp
                        expContainer.temp_cooldown = obj.FLIR.tempSeq;
                        obj.sc1Data.expList = [obj.sc1Data.expList expContainer];
                        
                        
                    end        %next speed     
                    
                    obj.alert.sound("start")
                    %Stop process
                    delete(poolObj);
                    obj.expAdd = 0;
                    obj.FLIR = thermalCamera(thermal_PORT);
                    obj.alert.sound("taskend");
                    fprintf("Scenario 1 finished\n")
 
                    
                case 2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SCENARIO 2
                    
                    %To complete 
                    
                case 3   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SCENARIO 3
  
                    %To complete
                    
                otherwise
                    fprintf(obj.alert.msg_scenario_error);
            end
            
        end
        
        function delete(obj)
           delete(obj.ADC);
           delete(obj.FLIR);
           delete(obj.controller);
           delete(obj.mCAP);
        end
    end
end

