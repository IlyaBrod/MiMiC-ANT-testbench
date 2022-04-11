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
classdef expData < handle
    %UNTITLED10 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name={}; %Experiment name
        id=0; %Experiment session id
        speed;%Set speed in %
        force;%Measured forces [time, Z, Y]
        temp; %Measured temperature [time, temperature]
        temp_cooldown; %Temperature of cooldown [time, temperature]
        voltage; %Measured voltage [time, voltage]
        current; %Measured current [time, current]
        legMotion; %[Time, Leg tip XYZ , ... , Leg base XYZ]
        groundMotion; %[Time, XYZ point Y increase & X increase order]
        date; %"23 Fev. 2022"
    end
    
    methods
        function obj = expData(name)
            %expData Construct an instance of this class
            obj.name = {name};
        end
    end
end

