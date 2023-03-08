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
classdef scData  < handle
    %UNTITLED11 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name={};
        id_count=0;
        mass;
        tempList; %[exp1temp, exp2temp...]
        loadList; %[exp1load, exp1load...]
        speedList; %[speed1,speed2]
        timeList; %[expTime1, expTime2] %length of the experimentation
        expList; %[exp1Container, exp2Container...]
    end
    
    methods
        function obj = scData(name)
            %scData Construct an instance of this class
            obj.name = {name};
        end
        
        function removeLast(obj)
            %remove Remove last experiment
            id = obj.expList(end).id;
            for exp=length(obj.expList):-1:1
                expID = obj.expList(exp).id;
                if(expID~=id)
                   obj.expList= obj.expList(1:exp+1);
                end
            end
            obj.id_count = obj.id_count -1;
            
        end
    end
end

