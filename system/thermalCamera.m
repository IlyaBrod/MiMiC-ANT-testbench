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

classdef thermalCamera < handle 
    %thermalCamera Handle FLIR thermal camera module
    %  Groups capture and configuration methods, run(dispInfo) for more
    %  possibilities.
    
    properties
        camera;
        frame;
        tempSeq;
    end
    
    properties (Access = private)
        portID;
    end
    
    methods
        function obj = thermalCamera(portID)
            %thermalCamera Construct thermal camera object
            %   PORT : serial port number (1,2,3,...)
            obj.portID = portID;
            obj.camera = videoinput('winvideo', obj.portID, 'Y16 _80x60');
            triggerconfig(obj.camera, 'manual'); %immediate
            %obj.camera.FramesPerTrigger = 1;
            %start(obj.camera);
        end
        
        function resetSequence(obj)
            obj.tempSeq = [];
        end
        
        function start(obj)
           start(obj.camera); 
        end
        
        function stop(obj)
           stop(obj.camera) 
        end
        
        function dispInfos(obj)
            cameraInfo = imaqhwinfo('winvideo', obj.portID);
            cameraInfo.SupportedFormats
            triggerinfo(obj.camera)  
        end
        
        function capture(obj)
            obj.frame = getsnapshot(obj.camera);
            obj.frame = double((obj.frame - 27315)) ./ 100;
        end
        
        function captureToBuffer(obj)
            shot = getsnapshot(obj.camera);
            shot = double((obj.frame - 27315)) ./ 100;
            obj.tempSeq = [obj.tempSeq, max(max(shot))];
            
        end
        function view(obj)
            figure;
            image(obj.frame,'CDataMapping','scaled')
            title('Thermal snaptshot')
            %set(gca,'ColorScale','log')
            colorbar
        end
        
        function lockUntil_max(obj,goalTemp,res)
            %lockUntil Lock the programm until a given temperature [°C]
            %res indicated the accepted error in °C.
            SETTINGS
            obj.resetSequence();
            currentTemp = goalTemp+1;
            tic;
            while(currentTemp>goalTemp+res || currentTemp<goalTemp-res)
                obj.capture;
                currentTemp = max(max(obj.frame));
                obj.tempSeq = [obj.tempSeq ; toc currentTemp];
                fprintf("Current Temp:%3.2f°C\n",currentTemp);
                pause(THREAD2_TECH);
            end
        end
        
        
        function lockUntil_mean(obj,goalTemp,res)
            %lockUntil Lock the programm until a given temperature [°C]
            %res indicated the accepted error in °C.
            currentTemp = goalTemp+1;
            while(currentTemp>goalTemp+res || currentTemp<goalTemp-res)
                obj.capture;
                currentTemp = mean(mean(obj.frame));
                fprintf("Current Temp:%3.2f°C\n",currentTemp);
            end
        end
        
        function reset(obj)
            stop(obj.camera);
            delete(obj.camera)
            obj.camera = videoinput('winvideo', obj.portID, 'Y16 _80x60');
            triggerconfig(obj.camera, 'manual'); %immediate
            start(obj.camera);
        end
        
        function preview(obj)
           %fastPreview Display capture witout temperature data
           delete(obj.camera)
           obj.camera = videoinput('winvideo', obj.portID, 'RGB24_80x60');
           triggerconfig(obj.camera, 'immediate');
           obj.camera.FramesPerTrigger = 1;
           
           figure;
           vidRes = obj.camera.VideoResolution;
           nBands = obj.camera.NumberOfBands;
           hImage = image( zeros(vidRes(2), vidRes(1), nBands) );
           preview(obj.camera,hImage)
            
           disp('Run a RESET after preview')
        end        

        function delete(obj)
           %delete Destructor method
           stop(obj.camera);
           delete(obj.camera);
        end
    end
end

