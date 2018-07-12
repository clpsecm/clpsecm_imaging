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
    start_location % scalar; start position (mm)
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

    function set_scan_distance(obj,start_location, length, resolution)
    % SET_SCAN_DISTANCE(start_location, length, resolution) in (mm)
        obj.start_location = start_location;
        obj.scanlength = length; 
        obj.resolution = resolution;
    end

    function lines = get_clpsecm_data(varargin)
    % lines = obj.GET_CLPSECM_DATA() Return SecmLines object from data.
        obj = varargin{1};
        SHEET_ONE = 2;
        DATA_OFFSET = 1;
        TICKS_COL = 'B';
        CURRENT_COL = 'C';

        %---data file path---%
        folderpwd = '../clpsecm_data/';
        filename = ['clpsecm_',obj.date,...
                    '_S', num2str(obj.sample_number,'%03d'),...
                    '_L', num2str(obj.nlines,'%02d'),...
                    '_',  num2str(obj.version,'%02d'),...
                    '.xlsx'];
        fpath = [folderpwd,filename];

        %---Find the scanning distance---%
        m = (obj.scanlength)     / (obj.resolution);
        s = (obj.start_location) / (obj.resolution);
        off = DATA_OFFSET;
        B   = TICKS_COL;
        C   = CURRENT_COL;

        ticks = xlsread( fpath, SHEET_ONE , ...
                    [B,num2str(s+off),':',B,num2str(s+m+off)] ); 
        ticks = ticks - (ticks(1)+ticks(end))/2;

        %---Generate line data---%
        currents = zeros(length(ticks), obj.nlines);          
        for L = 1:obj.nlines
            currents(:,L) = ... 
                xlsread( fpath, SHEET_ONE+(L-1) , ...
                    [C,num2str(s+off),':',C,num2str(s+m+off)] );  
        end
        
        if nargin == 2
            obj.params = varargin{2};
        else % read clpconfig
            obj.params = ProbeParams('config',ticks);
        end
        lines  = ScanLines(ticks, currents, obj.params);   
    end
end

end
