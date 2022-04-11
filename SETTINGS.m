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
%% Settings

CACHE_EXPORT = './Cache';

%Multicore data measure, needs at least 3 threads available
%if not, modify threading functions in MiMicAntBench.m  / run
%SAMPLING RATES
THREAD1_TECH = 0.004; % mocap camera 300Hz (majorized)
THREAD2_TECH = 1; %Thermal camera freq=8.6Hz (majorized)
%% Thermal camera
thermal_PORT = 1;
ROOM_TEMP = 25;
ROOM_TEMP_ERR = 0.5;
%% CONTROLLER
motor_VOLTAGE = 11.92;
arduino_PORT = "COM5";
arduino_BAUD = 9600;

%% NATIONAL INSTRUMENT settings
ni_RECORD_TIME = 30; %seconds (default time / auto modified)
ni_RECORD_RATE = 50000; %scans/s
ni_current_SENS = 0.1; %V/Amp  (Selected 100mV/A)
ni_current_PORT = 'ai2';
ni_sensZ_PORT = 'ai10';
%ni_voltage_PORT = 'ai3'; %need a step down to 10V
%% MOTION CAPTURE settings
% Custom marker name | Marker ID (given by gracking system).
% Nagative numbers for unlabelled markers
% Treadmill: TD (from left to right, from closer to base to farest one)
% Tibia:     TB (from down to top)
% Femur:     FE 
% Coxa:      CO 
% Torax:     TR 
% Calibrate  MK_CALIB (force plate calibration marker 1 on the top of each weight
% plus 3 markers on the corner, to be removed after.
% Order given by "labeled" menu in Qualisys software
mCAP_keySet = {'TD_1','TD_2','TD_3','TB_1', 'TB_2', 'TB_3', 'MK_CALIB_1', 'MK_CALIB_2', 'MK_CALIB_3', 'MK_CALIB_4'};
mCAP_valueSet = [-1 -2 -3 -1 5 6 1 -1 -2 -3];
mCAP_ID = containers.Map(mCAP_keySet,mCAP_valueSet); %Segments and point number

%% FORCE PLATE
%Calibration of Z SENSOR is done by run(-1) procedure. Should be
%executed everytime new motion capture calibration is done.
%Calibration of Y spring sensor, is done using a dynanometer and the
%procedure run(-2). The calibrated value can be stored in fPlate_STIFFNESS
%to avoid recalibration.
fPlate_calib_Xresolution = 5; %Number of subsections for the force plate calibration
fPlate_calib_Yresolution = 10;
fPlate_calib_minimal_layers = 4; %Minimum number of weight for calibration
                                 % Plus 0g weight
fPlate_STIFFNESS= 141.3953;%N/m (negative when not calibrated)

%Calibration structure, generated automatically.
%Possible to save and use the same several times, saving the
%MiMicAntBench.fPlate calibration
CALIB_VAR = 0.005; %Calibration variance
CALIB_DELAY = 0.5; %Calibration state change delay
fPlate_calib = struct('B',zeros(10,10),  ...
                    'xSize', .0,           ...
                    'ySize',.0,            ...
                    'xOffset',.0,           ...
                    'yOffset', .0,          ...
                    'Xresolution', fPlate_calib_Xresolution, ...
                    'Yresolution', fPlate_calib_Yresolution,...
                    'resetCalibration', true);

%% Leg parameters
LEG_MASS = 4.42/9.81*1e3; %mass in gramms
LEG_LOAD = 0;  %in gramms

%% SCENARIO 1: WALKING free treadmill
SC1_TIME = 5; %Total walking time in minutes
SC1_SPEED_STEPS = 2; %Walking speed steps (put 2 if only one speed available + static)
SC1_MAX_SPEED = 30; % 0-100% (ALWAYS > MIN SPEED)
SC1_MIN_SPEED = 0; % 0-100% -> IF ONLY 1 SPEED AVAILABLE SET HERE
