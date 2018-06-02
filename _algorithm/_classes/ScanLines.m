classdef ScanLines
    properties
        nlines;     %-scalar; number of lines
        nmeasures;  %-scalar; number of measurements in each line
        currents;    %-(nmeasures,nlines); current data
        angles;     %-(1, nlines);  angle of scan lines
        xticks      %-(nmeasures, 1); xticks of measuemrnt (mm)
    end
    
    methods 
        function obj = ScanLines(angles)
            obj.nlines = length(angles);
            obj.angles = angles;
        end
    end
end