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

classdef resultTable < handle
    %UNTITLED10 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        data;
    end
    
    methods
        function obj = resultTable()
            %resultTable Construct result table and init

            varTypes = {'string','string','string','string'};
            sz = [12 4];
            
            resultTable_rows = {'Mass [kg]' 'Max|Mean|Min. Power [W]' 'Max. speed [cm/s]' ...
                    'Stride length [mm]' 'Stride max. freq. [Hz]' 'Innefective stance [%]' ...
                    'Step accuracy [mm]' 'Max tangent stress [N]' 'Max. normal stress [N]' ...
                    'CoT' 'Temp. inc. [°/min]' 'Comment'};

                resultTable_columns = {'Leg' 'Hexapod' 'Scale_factor' 'Reference'};

                obj.data = table('Size',sz,'VariableTypes',varTypes); % Table init
                obj.data.Row = resultTable_rows;
                obj.data.Properties.VariableNames = resultTable_columns;
        end
        
        function get(obj,x,y)
            %get Return an value on the ligne x, column y
            disp(obj.data(x,y));
        end
        
        function obj=set(obj,x,y,value)
            %set Change an value on the ligne x, column y to value
            obj.data(x,y) = {value};
        end
        
        function importCSV(obj,filename)
            %importCSV Give a filename to import precalculated table data
            %  CSV file format: separation ',', no value '-'
            idata = readtable(filename,'ReadRowNames',true);
            for x=1:size(obj.data)*[1 0]' 
                for y=1:size(obj.data)*[0 1]'
                   if( string(table2cell(idata(x,y))) ~= "")
                       obj.data(x,y) = idata(x,y);
                   end
                end
            end
        disp("Importe done");    
        end
        
        
        function export(obj,filename)
           %export Save the table data to a csv file
           writetable(obj.data,strcat(filename,'.csv'),'Delimiter',',','QuoteStrings',true,'WriteRowNames',true)
           disp(strcat('Result file: ',filename,'.csv saved'))
        end
        
    end
end



