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
classdef alertBox < soundPlayer
        %alertBox Manage alerts sounds and text

    properties
        
        %TODO : organize better all this :)
        %TEXT MESSAGES                
        msg_scenario_error = join(["Incorrect scenario selected\n"
                            "Please refer to the documentation\n"]);
        %SCENARIO 1 MESSAGES
        msg_scenario_1_warn = join(["SCENARIO 1 selected\n"
                   "Please, process the following checklist verification:\n"
                   "\t1- Spring detached from the calibrated treadmill\n"
                   "\t2- Leg vertical axe unlocked with stop at the lowest position\n"
                   "\t3- 3D marker placed (leg tip, base), system calibrated\n"
                   "\t4- Check the thermal camera position\n"
                   "\t5- Put the current sensor ON at 100mV/A\n"]);
       
        
        %CALIBRATION MESSAGES (SCENARIO -1)
        msg_calibration_warn = join(["Starting force plate calibration procedure\n"
                    "\t1- At least 3 weight sould be used to calibrate the system"
                    "with a value covering the force range needed\n"
                    "\t2- Each weight should be equipped of a 3D marker on the center\n"
                    "\t3- Be sure to identificate the 3D marker ID in SETTINGS (MK_CALIB)\n"
                    "PROCEDURE- The masses are going to be placed one by one on the treadmill"
                    "and moved when idicated, in order to cover all the plate surface\n"]);
      msg_calibration_init = "Remove all weight from the treadmill\n";
      msg_calibration_corner1 = "Place markers 2, 3, 4 on the plate corners\n";
      msg_calibration_corner2 = "Remove all markers from the corners\n";
      msg_calibration_addmass = "Please take a new weight value\n";
      msg_calibration_number = "Enter available weight number: ";
      msg_calibration_process = join(["\tINSTRUCTIONS: Add enough positions to complete the pattern\n"
                                    "Witout lifting the weight\n"]);
      msg_calibration_addposition = join(["Move the weight to an"
                    " other position on the treadmill\n"]);
      msg_calibration_value = 'Indicate mass value [g]: ';
      msg_calibration_zero = "Keep clean the plate during 0g calibration\n";
      
      %CALIBRATION Y MESSAGES (SCENARIO -2)
      msg_calibration2_warn =  join(["Starting force plate Y axe calibration procedure\n"
                    "\t1- You should be equipped with a dynamometer on a sliding card"
                    "on the level of the threadmill top surface\n"
                    "\t2- Constrain the threadmill with the Y spring using the attach point\n"
                    "\t3- Use the top surface attach part, to join the dynamometer and the mat\n"
                    "\t4- The threadmill should be equipped by TD markers (see SETTINGS)\n"
                    "PROCEDURE- You have to pull dynamometer and accordingly the mat."
                    "Then specify the dynamometer value.\n"]);
      
    end
    
    methods
        function obj = alertBox(~)
            
        end
    end
end

