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

classdef robotController
    %RobotController manage serial connection to start code executing.
    % Set MANUAL=1 if the speed should be switched manually on demand
    % Serial data format: 
    % Start : 'SE' + char(scenario number)
    % Stop : 'SE' + char(0)
    % Speed in 0-100% : 'SP' + char(speed value (0-100))
    properties
      COM="COM1";
      MANUAL = 0;
      BAUDRATE=9600;
      myController;
    end
    
    methods
        function obj = robotController(COM,BAUDRATE,MANUAL)
           obj.MANUAL = MANUAL;
           obj.COM = COM;
           obj.BAUDRATE = BAUDRATE;
           obj.myController = serial(obj.COM,'BaudRate',obj.BAUDRATE,'terminator','CR');
           fopen(obj.myController);
        end

        function serialReset(obj)
            fclose(obj.myController);
            fopen(obj.myController);
        end
        
        function serialStop(obj)
           fclose(obj.myController); 
        end
        
        function start(obj,scenario)
           %fwrite(obj.myController,'end');
           fwrite(obj.myController,['SE' char(scenario)]);
           
        end
           
        function stop(obj)
           fwrite(obj.myController,['SE' char(0)]);
        end
        
        function setSpeed(obj,ratio)
           %setSpeed configure controller speed, if possible, using a
           %ratio from 0 to 100 (min-max).
           if(obj.MANUAL==1)
               fprintf("Set speed MANUALLY to %d %% \n",ratio);
               wait;
           else
               fwrite(obj.myController,['SP' char(ratio)]);
           end
           
        end
        
        function delete(obj)
             fclose(obj.myController); 
        end
     end

end

