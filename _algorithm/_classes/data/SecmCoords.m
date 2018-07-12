classdef SecmCoords < handle & matlab.mixin.Copyable
% SECMCOORDS Creates the coordinate system of clp-secm imaging system
%   obj = SECMCOORDS(ticks) Creates coordinates, ticks should be
%   equi-stepsized sequence. The algorithm automatically centering 
%   the ticks.
%
%   The class SECMCOORDS and its childs are handle objects.
%
% SECMCOORDS methods:
%   DOWNSAMPLE - Downsample the resolution of coordinates.
%   FFT        - Fourier transform in SECM coordinates.
%   IFFT       - Inverse Fourier transform in SECM coordinates.
%
% SECMCOORDS properties:
%   NMEASURES  - Number of pixls in x(y) directions.
%   TICKS      - The ticks in both xy direction in (mm)
%   RESOLUTION - Resolution of coordinates system in (mm)
%   XYLIM      - The boundary of coordinate in (mm) 
% 
% See also SECMIMAGE, SCANLINES

properties (SetAccess = private)
    xylim      % vector(2,1); two ends of scan on both direction.
    ticks      % (nmeasures, 1); The ticks in both xy direction.
    resolution % scalar; resolution of coordinates in (mm)
    nmeasures  % scalar; number of ticks in x(y) directions.
end


properties (Access = protected)
    Ceff  % [nmeasures,nmeasures]; Effective scannning area
end

methods
    function obj = SecmCoords(ticks)
        % obj = SECMCOORDS(ticks) Creates coordinates, ticks should 
        % be an equi-stepsized sequence. The algorithm automatically
        % centering the ticks.
        if isrow(ticks); ticks = ticks'; end
        obj.resolution = mean(ticks(2:end)-ticks(1:end-1));
        obj.ticks = ticks - (ticks(1)+ticks(end))/2;
        obj.set_coords_Ceff();
    end

    function downsample(obj,downrate)
        % obj.DOWNSAMPLE(downrate) Downsample the coordinate by an 
        % integer downsample rate.
        if mod(downrate,1)~=0
            error('Input downrate should be integer'); 
        end

        idxlist = 1:downrate:length(obj.ticks);
        if obj.ticks(end) ~= obj.ticks(idxlist(end))
            error('Downsampling throws away end tail of data.')
        end
        
        obj.resolution = obj.resolution*downrate;
        obj.ticks = obj.ticks(idxlist);
        obj.set_coords_Ceff();
    end
end

methods (Static)
    function fdata = fft(data)
        % fdata = obj.FFT(data) Static; Fourier transform.
        fdata = fft2(ifftshift(data));
    end
    function data = ifft(fdata)
        % data = obj.IFFT(fdata) Static; Inverse Fourier transform.
        data = fftshift(ifft2(fdata));
    end
end   

methods (Access = private)
    function set_coords_Ceff(obj)
        % Set coordinate dimension and Ceff
        res = obj.resolution;
        obj.nmeasures = length(obj.ticks) ;
        obj.xylim = [min(obj.ticks)-res/2, max(obj.ticks)+res/2];

        N = obj.nmeasures;
        T2 = (obj.ticks).^2;
        R = max(abs(obj.xylim));
        obj.Ceff = repmat(T2,[1,N]) + repmat(T2',[N,1]) <= R^2;
    end
end
    
end % classdef