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
%% Test bench running programm :
    % Init sensors
    % Run the tests and display results (depending the scenario)
%% Scenario selection
% 1- Walk 
% (Stride length ; Stride max freq ; Innefective stance ;
% step accuracy ; CoT ; Temp inc.)
% 2- Walk with resistance
% (Max. power ; Max tan ; Temp inc. ;
% 3- Weight lifting 
% (Max norm stress ; Temp. inc.)

SCENARIO = 1;


%% Init process
bench = MiMicAntBench();

% Helps to identify system points
%bench.mCAP.listPoints 

% Check Thermal camera position
%bench.FLIR.preview
%bench.FLIR.reset

%Run scenario
bench.run(SCENARIO);

%Remove experiment session
bench.sc1Data.removeLast

%Add to previous experiment session
bench.addExperiment;
bench.run(1);
%bench.cancelAddExperiment()


%% Ending experiment
%delete(bench)



