classdef ScanLines < SecmCoords
% SCANLINES Creates or modify the current data of CLP-SECM lines.
%
% obj = SCANLINES(ticks,currents,params) Creates scan line object with
% sepcified ticks, currents and params with format:
%       ticks:    vector(nmeasures,1); CLP scans distance labels.
%       currents: matrix(nmeasures,nlines); CLP current measures.
%       params:   object CalibParams; CLP probe properties.
%
% obj = SCANLINES(ticks,currents) Create scan line object using parametric
% setting in clpconfig.m
%
% obj = SCANLINES(obj1) Copy construtor.
%
% SCANLINES is a handle object.
%
% SCANLINES methods:
%   PLUS   -  Overload '+'. return new SCANLINES w/ summed currents.
%   MINUS  -  Overload '-', return new SCANLINES w/ subtracted currents.
%   NORM   -  Return Frobenius norm of currents.
%   INPROD -  Return inner product of currents.
%   CONV   -  Convolve currents w/ input kernel.
%   MULT   -  Pointwise multiply currents w/ input multiplier(s).
%   SHIFT  -  Shift currents data by input distance in (mm).
%   BACK_PROJECT - Backproject of currents to SECM image base on params
%   DOWNSAMPLE   - Decrease resolution of current data by integer rate.
%   PLOT_LINES   - Plot current data (with selected angles in deg).
%  
% SCANLINES public fields:
%   NLINES   -  scalar; number of lines.
%   CURRENTS -  vector(nmeasure,nlines); current data. 
%   PARAMS   -  object ProbeParams; CLP probe properties.
%
% See also SECMIMAGE, SECMCOORDS, PROBEPARAMS, CLPCONFIG

properties 
    nlines    % scalar; number of lines.
    currents  % (nmeasures,nlines); current data.
    params    % object ProbeParams; CLP probe properties.
end

methods 
    function obj = ScanLines(varargin)
    % Creates scan line object with sepcified angles, distance and current data. 
        if nargin == 1 % Copy Contstructor
            obj1 = varargin{1};
            ticks = obj1.ticks;
            nlines = obj1.nlines;
            params = obj1.params;
            currents = obj1.currents;
        elseif nargin == 2 % Input type (ticks,currents)
            ticks = varargin{1};
            currents = varargin{2};
            nlines = size(currents,2);
            params = ProbeParams('config',ticks);     
        elseif nargin == 3 % Input type (ticks,currents,params)
            ticks = varargin{1}; 
            currents = varargin{2};
            nlines = size(currents,2);
            params = varargin{3}; 
        end
        obj = obj@SecmCoords(ticks);
        obj.nlines = nlines;
        obj.params = params;
        obj.currents = currents;
    end

    function obj = plus(obj1,obj2)
    % obj = obj1 + obj2 Return new SCANLINES w/ summed currents.  
        check_issame(obj1,obj2)
        obj=ScanLines(obj1.ticks, obj1.currents+obj2.currents, obj1.params);             
    end

    function obj = minus(obj1,obj2)
    % obj = obj1 - obj2 Return new SCANLINES w/ subtracted currents.
        check_issame(obj1,obj2)
        obj=ScanLines(obj1.ticks, obj1.currents-obj2.currents, obj1.params);              
    end

    function n = norm(obj); n = norm(obj.currents,'fro'); end
    % n = obj.NORMS() Return frobenius norm of currents

    function n = inprod(obj1,obj2) 
    % n = obj.INPROD() Return frobenius norm of currents
        check_issame(obj1,obj2);
        n = sum(sum(obj1.currents.*obj2.currents)); 
    end
        
    function conv(obj,psf)
    % obj.CONV(psf) Convolute currents w/kernels. 
        if isscalar(psf)
            obj.currents = obj.currents * psf;
        else % kernel is vector of length nmeasures
            F  = @(line)  SecmCoords.fft(line);
            iF = @(fline) SecmCoords.ifft(fline);
            fk = F(psf);
            for I = 1:obj.nlines
                obj.currents(:,I) = real(iF(F(obj.currents(:,I)).*fk));
            end
        end
    end

    function mult(obj,multiplier)
    % obj.MULT(multiplier) Point-product currents w/ multiplier. 
        obj.currents = obj.currents .* multiplier;
    end

    function shift(obj,s)
    % obj.SHIFT(s) Shift i-th currents by s(i) for all currents.
        d = s/obj.resolution;
        for I = 1:obj.nlines
            if mod(d(I),1) == 0 % d is integer
                obj.currents(:,I) = shift(obj.currents(:,I),d);
            else % d is not integer
                dl = floor(d(I));
                dh = ceil(d(I));
                obj.currents(:,I) = ...
                    (dh-d(I)) * shift(obj.currents(:,I),dl) + ...
                    (d(I)-dl) * shift(obj.currents(:,I),dh) ;
            end
        end
    end

    function secmImage = back_project(obj)
    % secmImage = obj.BACK_PROJECT() Back project currents to image.
    % BACK_PROJECT is the adjoint operator of line projection. 
    %
    % See also SECMIMAGE.LINE_PROJECT
        N = obj.nmeasures;
        obj1 = ScanLines(obj);
        
        % Reverse (adjoint) convolution with kernel 
        if obj1.params.psf.isactive
            p = obj1.params.psf.func(obj1.params.psf.value);
            obj1.conv(flip(p));
            if norm(p) == 0
                error('Point spread function should be nonzero');
            end
        end
        
        % Multiply currents data by multiplier
        if obj1.params.intensity.isactive
            m = obj1.params.intensity.func(obj1.params.intensity.value);
            obj1.mult(m);
        end

        % Reverse (adjoint) shift all lines
        if obj1.params.shifts.isactive
            s = obj1.params.shifts.func(obj1.params.shifts.value);
            obj1.shift(-s);
        end

        % Back project using the scanning angles
        angles = obj1.params.angles.func(obj1.params.angles.value);
        img = zeros(N,N);
        for I = 1:obj1.nlines
            img = img + SecmImage.rotate_image(...
                repmat(flip(obj1.currents(:,I)),[1,N]), angles(I));
        end
        secmImage = SecmImage(obj1.ticks,img);
        secmImage.image = (secmImage.image).*(secmImage.Ceff);
    end

    function downsample(varargin)
    % obj.DOWNSAMPLE(downrate) Downsample lines by interger rate.
    % obj.DOWNSAMPLE(downrate,params) Downsample and assign parameter.
        obj = varargin{1};
        downrate = varargin{2};
        idxlist = 1:downrate:length(obj.ticks);
        obj.currents = obj.currents(idxlist,:);
        downsample@SecmCoords(obj,downrate);
        if nargin == 3
            obj.params = varargin{3};
        else
            obj.params = ProbeParams('config',obj.ticks);
        end
    end

    function plot_lines(obj,angles_subset) 
    % obj.PLOT_LINES() Plot all current data
    % obj.PLOT_LINES(angles_subset) Plot lines with input angles.
        angles = obj.params.angles.value;
        if nargin == 1; angles_subset = angles; end

        hold on; box on; legend_info = {};
        if iscolumn(angles_subset); angles_subset = angles_subset'; end
        for angle = angles_subset
            idx = find(angles == angle);
            if idx
                plot(obj.ticks, obj.currents(:,idx));
                legend_info{end+1} = ['scan angle: ',num2str(angle,'%.2f'),'º'];
            end
        end
        legend(legend_info,'location','northwest'); 
        xlabel('Distance/mm'); xtickformat('%d');
        ylabel('Current/A'); ytickformat('%1.2f');
    end
end
    
methods (Access = private)
    function check_issame(obj1,obj2)
        if ~isequal(obj1.ticks,obj2.ticks)
            error('SecmImage convolution available only in same coordinates')
        end
    end
end

end