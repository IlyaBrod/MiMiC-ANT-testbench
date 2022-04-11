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

classdef soundPlayer < handle
    %soundPlayer Module to play sounds using VLC
    %Warning: Add VLC to system path
    
    properties
        state=0;
    end
    
    methods
        function obj = soundPlayer(~)
            obj.state=0;
        end
        
        function sound(obj,ID)
            %Play alert sound depending the provided alert ID
            % 1- camera
            % 2- alert
            % 3- bip
            % 4- counter
            % 5- taskend
            % 6- message
            % 7- start
            % 8- stop
            switch (ID)
                case "camera"
                    !vlc ./Sounds/s_camera.oga --qt-start-minimized --play-and-exit&
                case "alert"
                    !vlc ./Sounds/s_alert.oga --qt-start-minimized --play-and-exit&
                case "bip"
                    !vlc ./Sounds/s_bip.wav --qt-start-minimized --play-and-exit&
                case "counter"
                    !vlc ./Sounds/s_counter.oga --qt-start-minimized --play-and-exit&
                case "taskend"
                    !vlc ./Sounds/s_taskend.wav --qt-start-minimized --play-and-exit&
                case "message"
                    !vlc ./Sounds/s_message.oga --qt-start-minimized --play-and-exit&
                case "start"
                    !vlc ./Sounds/s_start.oga --qt-start-minimized --play-and-exit&
                case "stop"
                    !vlc ./Sounds/s_stop.oga --qt-start-minimized --play-and-exit&
            end
        end
    end
end

