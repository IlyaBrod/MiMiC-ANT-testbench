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

classdef MotionCapture < handle 
    %motionCapture Manage Qualisys motion capture system
    
    properties
        ip;
        labels;
    end
    
    properties (Access = private)
        data;
        dataNoLabel;
        frameinfo;
        
    end

    methods (Access = private)
        function getData(obj)
            [obj.frameinfo,obj.data, obj.dataNoLabel] = QCM;
             %Get rid of the inital state
             if(obj.frameinfo(2) == 0)
                 [obj.frameinfo,obj.data] = QCM;
             end
        end
        
    end
    
    methods
        function obj = MotionCapture(~)
            %MotionCapture Initialize motion capture parameters
            % to 3d Euler
            
            obj.ip = '127.0.0.1';
            QCM('connect', obj.ip, 'frameinfo', '3D 3DNoLabels'); %3DNoLabels
            
            %Groups init
            labels_int = QCM('6doflabels');
            obj.labels=deblank(string(labels_int)); % Convert to string array
            
            %Display init

        end
       
        function updateLabels(obj)
            %updateLabels update the local 3D objects name list
            labels_int = QCM('6doflabels');
            obj.labels=deblank(string(labels_int)); % Convert to string array
        end
        
        function listPoints(obj)
            obj.getData;
            fprintf("Points with label\n")
            fprintf("[X\tY\tZ\tResidual]\n")
            obj.data
            fprintf("Points w/o label:\n")
            fprintf("[X\tY\tZ\tResidual]\n")
            obj.dataNoLabel
        end
        
        function [x,y,z] = getPoint(obj,pointID)
           %getPoint return the 3d coordinates of the asked labeled point
           % check label container to know the possible values (SETTINGS)    
           
           x = NaN;
           y = NaN;
           z = NaN;
           
           try
               obj.getData;
               
               if(pointID>=0)
                selected_data = obj.data;
               else
                selected_data = obj.dataNoLabel;
               end

               x = selected_data(abs(pointID),1);
               y = selected_data(abs(pointID),2);
               z = selected_data(abs(pointID),3);
           end
               
        end
        
        function delete(~)
            QCM('disconnect');
            clear mex
        end
    end
end

