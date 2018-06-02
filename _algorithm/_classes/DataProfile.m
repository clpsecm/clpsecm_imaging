classdef DataProfile
    properties
        date = '053018';     %-string; MMDDYY
        sample_number = 0;   %-scalar; sample number
        nlines = 4;          %-scalar; number of scan lines
        version = 1;         %-scalar; versions of same sample,nlines on same date
        scan_start_location = 100; %-scalar; start position (mm)
        scan_length = 700;         %-scalar; length of scan in (mm)
        scan_resolution = 0.1;     %-scalar; resoluiton of scan in (mm)
        angles = [0,40,80,120];    %-(1,nlines); angles of scan lines
    end
end
