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

classdef NationalInstrument < handle
    %NationalInstrument class for manage DAC device
    %   Add channels, record, all grouped !
    
    properties
        port;
        DAC;
        channels= [];
        channel_names = [];
        data;
        time;
        threadIT=0;
    end

    properties (Access = private)
       ls; 
    end
    methods
        function obj = NationalInstrument(ni_RECORD_RATE,ni_RECORD_TIME)
            devices = daq.getDevices;
            obj.port = devices(1).ID;
            obj.DAC = daq.createSession('ni');
            obj.DAC.Rate = ni_RECORD_RATE;
            obj.DAC.DurationInSeconds = ni_RECORD_TIME;
        end
        
        function obj = addChannel(obj,customName,portName)
            newChannel = DACinput(customName,obj.DAC,obj.port,portName,'Voltage');
            obj.channels = [obj.channels; newChannel];
            obj.channel_names = [obj.channel_names; string(customName)];
        end
        
        function id = getIndex(obj,customName)
            %getIndex Return the channel index name for data extraction
            id = find(ismember(obj.channel_names,string(customName)),1);
        end
        
        function output = read(obj)
           %read Return instant value acquired by continues recording
           [idata, ~] = obj.DAC.inputSingleScan;
           output = idata;
        end
        
        function output = record(obj)
            %record Start capture for the configured settings
            %return [time, data]
            [obj.data,obj.time] = obj.DAC.startForeground;
            output = [obj.time,obj.data];
        end
        
        function output = recordSeconds(obj,time)
            %recordSeconds Start capture for a precise time lapse in sec.
            %return [time, data]
            saveTime = obj.DAC.DurationInSeconds;
            obj.DAC.DurationInSeconds = time;
            [obj.data,obj.time] = obj.DAC.startForeground;
            obj.DAC.DurationInSeconds = saveTime;
            output = [obj.time,obj.data];
        end
         
        function recordThread(obj,time)
            %record Start capture in background
            %Data saved locally for threading
            %return None
            prepare(obj.DAC);
            obj.threadIT=0;
            obj.data = zeros(obj.DAC.Rate*time,length(obj.channel_names));
            obj.time = zeros(obj.DAC.Rate*time,1);
            %func = @NationalInstrument.plotData;
            obj.ls = addlistener(obj.DAC,'DataAvailable',@(scr,event)obj.thread(obj,scr,event,1)); 
            obj.DAC.startBackground;
            fprintf("\tADC thread started\n");
        end
        
        function stopThread(obj)
           delete(obj.ls);
        end
        
        function thread(obj,src,event,varargin)
            it = obj.threadIT;
            len = 5000;
            %size(varargin{1}.Data)
            %size(varargin{1}.TimeStamps)
            obj.data(it*len+1:(it+1)*len,:) = varargin{1}.Data;% [obj.data ; varargin{1}.Data];
            obj.time(it*len+1:(it+1)*len,:) = varargin{1}.TimeStamps;%[obj.time ; varargin{1}.TimeStamps];
            obj.threadIT=obj.threadIT+1;
            %plot(varargin{1}.TimeStamps,varargin{1}.Data)
        end
        
         function plot(obj)
            figure;
            
            plot(obj.time,obj.data)
            xlabel('Time (secs)');
            ylabel('Voltage');
            title('DAQ raw data')
            legendStr = cell(1,length(obj.channels));
            for i=1:numel(obj.channels)
                legendStr(1,i) = {obj.channels(i).name};
            end
         end
         
         function setRecordTime(obj,time)
             obj.DAC.DurationInSeconds = time;
         end
    end
    
end

