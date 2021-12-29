classdef DataSpec < handle
% DATASPEC Register profile of data file name and its spec.
%
% obj = DATASPEC(date, sample_number, nlines, version) Specify 
%       date           - string; MMDDYY
%       sample_number  - serial number of sample. 
%       nlines         - number of scan lines.
%       version        - version number.
%
% DATASPEC methods:
%   SET_SCAN_DISTANCE - Set scan start_location, length, resolution in (mm)
%   GET_CLPSECM_DATA  - Return line data object from file and clpconfig 
% 
% See also SCANLINES, CLPCONFIG

properties (SetAccess = private)
    date           % string; MMDDYY
    sample_number  % scalar; sample number
    nlines         % scalar; number of scan lines
    version        % scalar; version number
    start_time     % scalar; start time (s)
    start_location % scalar; start position (mm)
    nskip          % scalar; downsampling rate
    scanlength     % scalar; length of scan in (mm)
    resolution     % scalar; resoluiton of scan in (mm)
    params         % ProbeParam; parameters read from clpconfig
end

methods
    function obj = DataSpec(date, sample_number, nlines, version)
    % Construct object with chosen file
        obj.date = date;
        obj.sample_number = sample_number;
        obj.nlines = nlines;
        obj.version = version;
    end

    function set_scan_distance(obj, resolution, length, start_time, nskip)
    % SET_SCAN_DISTANCE(start_location, length, resolution) in (mm)
        obj.start_time = start_time;
        obj.start_location = start_time * resolution;
        obj.scanlength = length; 
        obj.resolution = resolution;
        obj.nskip = nskip;
    end

    function lines = get_clpsecm_data(varargin)
    % lines = obj.GET_CLPSECM_DATA() Return SecmLines object from data.
        obj = varargin{1};
        SHEET_OFFSET = 2;
        DATA_OFFSET = 1;
        TIME_COL = 'B';
        CURRENT_COL = 'C';
        TIME_PERIOD = 0.1; % Current recoding period in (sec)
        SAMPLE_PERIOD = 1; % Sample time period in (sec)
        
        %---data file path---%
        folderpwd = '../clpsecm_data/';
        filename = ['clpsecm_', obj.date,...
                    '_S', num2str(obj.sample_number,'%03d'),...
                    '_L', num2str(obj.nlines,'%02d'),...
                    '_',  num2str(obj.version,'%02d'),...
                    '.xlsx'];
        fpath = [folderpwd,filename];

        %---Find the scanning distance---%
        s_t  = SAMPLE_PERIOD / TIME_PERIOD; 
        m = floor(s_t * obj.scanlength / obj.resolution);
        s = s_t * obj.start_time;
        doff = DATA_OFFSET;
        B    = TIME_COL;
        C    = CURRENT_COL;


        % Generate ticks and line data
        sampletime = xlsread( fpath, SHEET_OFFSET, ...
                [B,num2str(s+doff),':',B,num2str(s+m+doff)] ); 
        currents = zeros(numel(sampletime), obj.nlines);          
        for L = 1:obj.nlines
            currents(:,L) = ... 
                xlsread( fpath, SHEET_OFFSET+(L-1) , ...
                    [C,num2str(s+doff),':',C,num2str(s+m+doff)] );
        end
        
        % Find sampled current and ticks, center and calculate the ticks
        sampletime = sampletime(1: obj.nskip :end);
        currents = currents(1: obj.nskip :end, :);
        
        % Remove tail zeros from insufficient scan length
        sampletime = sampletime(1 : find(sampletime,1,'last'));
        currents = currents(1:length(sampletime), :);
        
        ticks = (sampletime - obj.start_time) * obj.resolution;
        ticks = ticks - (ticks(1)+ticks(end))/2; 
        
        
        %---Generate ScanLines data---%
        if nargin == 2
            obj.params = varargin{2};
        else % read clpconfig
            obj.params = ProbeParams('config',ticks);
        end
        lines  = ScanLines(ticks, currents, obj.params);   
    end
end

end
