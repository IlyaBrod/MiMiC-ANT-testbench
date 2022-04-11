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
classdef forcePlate < handle
    %forcePlate Init force plate sensor (Normal and tangential sensors)
    
    properties
        calibration;
        stiffness;
    end
    
    properties (Access = private)
        xinterval;
        yinterval;
        calibrated;
    end
    
    methods
        function obj = forcePlate(~)
            %forcePlate Construct forcePlate sensor
            obj.stiffness = 0.5;
            obj.calibrated = false;
        end
        
        function init_calibrationZ(obj, calib)
            %Initialize calibration values
            obj.calibration = calib;
            obj.xinterval = calib.xSize / calib.Xresolution;
            obj.yinterval = calib.ySize / calib.Yresolution;
            if(calib.resetCalibration==false)
                obj.calibrated = true;
            end
        end
        
        
        function force = getForce_norm(obj,x,y,voltage)
            %assert(obj.calibrated==true,"Error: System not calibrated")
            %[i,j] = obj.getIndex(x,y);
            X = [1 voltage x y];
            force = (X*obj.calibration.B).*9.80665;
        end
        
        function force = getForce_tan(obj,distance)
            %getDorce_tan Return tangentiel force in Newton
            %distance is given in meters.
            force = distance*obj.stiffness;
            
        end
        
        function [i, j] = getIndex(obj,x,y)
            i = fix((x-obj.calibration.xOffset)/obj.xinterval)+1;
            j = fix((y-obj.calibration.yOffset)/obj.yinterval)+1;
        end
        
        function [x,y] = reverseIndex(obj, i,j)
            x = (i-1)*obj.xinterval+obj.calibration.xOffset;%bki * obj.xinterval -1; 
            y = (j-1)*obj.yinterval+obj.calibration.yOffset;;%bkj * obj.yinterval -1;
        end
        
    end
end

