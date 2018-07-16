classdef SecmImage < SecmCoords
% SECMIMAGE Creates or modify image with SECM coordinate system. The image
% will always supported only on Ceff.
% 
%   obj = SECMIMAGE(ticks) Creates object with assigned coordinate
%   and empty image.
% 
%   obj = SECMIMAGE(ticks,image) Creates object with with assigned 
%   coordinates and input image.
%
% SECMIMAGE is a handle object.
%   
% SECMIMAGE methods:
%   PLUS   -  Overload '+'; return new SECMIMAGE w/summed image.
%   MINUS  -  Overload '-'; return new SECMIMAGE w/subtracted image.          
%   MTIMES -  Overload '*'; return new SECMIMAGE w/convolved images.
%   SOFT   -  Return obj with soft-thresholded image by lda.
%   POS    -  Return obj with image w/o negative entries
%   SUM    -  Return sum of image values. 
%   NORM   -  Return Frobenious norm of image.
%   MAX    -  Return max value in image.
%   INPROD -  Return inner product between image and input mask.
%   ROTATE -  Rotate image by an input angle counterclockwise.
%   SET_IMAGE    - Assign input image.
%   LINE_PROJECT - Project image to lines of given slope angles.
%   DOWNSAMPLE   - Decrease resolution of image. 
%   DRAW_IMAGE   - Draw image with inherent coordinate system
%   
% SECMIMAGE public fields: 
%   IMAGE      - Discrete image with inherit coordinates. 
%   NMEASURES  - Number of pixls in x(y) directions.
%   TICKS      - The ticks in both xy direction in (mm)
%   RESOLUTION - Resolution of coordinates system in (mm)
%   XYLIM      - The boundary of coordinate in (mm) 
%
% See also SECMCOORDS, SCANLINES, DICTPROFILE, SPARSEMAP

properties
    image  % matrix(nmeasures,nmeasures); The image.
end

methods
    function obj = SecmImage(ticks,img)
        % Constuct object SECMIMAGE of input coordinates.
        obj = obj@SecmCoords(ticks); 
        N = obj.nmeasures;
        switch nargin
            case 1; obj.image = zeros(N,N);
            case 2; obj.image = (obj.Ceff).*img;
        end
    end

    function obj = plus(obj1,obj2)
    % obj = obj1 + obj2 Return a new SECMIMAGE with summed image 
        check_issame(obj1,obj2);
        obj = SecmImage(obj1.ticks, obj1.image + obj2.image);
    end

    function obj = minus(obj1,obj2)
    % obj = obj1 - obj2 Return a new SECMIMAGE with subtracted image
        check_issame(obj1,obj2);
        obj = SecmImage(obj1.ticks, obj1.image - obj2.image);
    end

    function obj = mtimes(obj1,obj2)
    % obj = obj1 * obj2 Return a new SECMIMAGE with convolution of images.
        if isa(obj1,'SecmImage')
            obj = SecmImage(obj1.ticks);
            if isa(obj2,'SecmImage')
                check_issame(obj1,obj2);
                fimage = SecmCoords.fft(obj1.image).* ...
                         SecmCoords.fft(obj2.image);
                obj.image = real(SecmCoords.ifft(fimage));
            else % obj2 is scalar
                obj.image = obj1.image * obj2;
            end
        else % obj1 is scalar
            obj = SecmImage(obj2.ticks);
            obj.image = obj2.image * obj1;
        end        
    end

    function obj = soft(obj,lda); obj.image = soft(obj.image,lda); end
    % obj = obj.SOFT(lda) Return soft-thresholded image w/ threshold lda.

    function obj = pos(obj); obj.image = max(obj.image,0); end
    % obj = obj.POS(); Return obj with image w/o negative entries.

    function s = sum(obj);  s = sum(sum(obj.image)); end
    % s = obj.SUM(); Return sum of image entries.

    function s = norm(obj); s = norm(obj.image,'fro'); end 
    % s = obj.NORM(); Return Frobenious norm of image.

    function s = max(obj); s = max(max(obj.image)); end
    % s = obj.MAX(); Return max value of image.

    function s = inprod(obj1,obj2) 
    % s = obj.INPROD(mask); Return inner product of two image
        check_issame(obj1,obj2);
        s=sum(sum(obj1.image.*obj2.image)); 
    end 
    

    function set_image(obj,img); obj.image = (obj.Ceff).*img; end
    % obj.SET_IMAGE(img); Set image as img

    function obj = rotate(obj,angle)
    % obj = obj.ROTATE(angles); Rotate image by angle counterclockwise
        obj.image = obj.Ceff.*(obj.rotate_image(obj.image,angle));
    end
    
    function downsample(obj,downrate)
    % obj.DOWNSAMPLE(downrate); Downsample image by integer rate. 
        idxlist = 1:downrate:length(obj.ticks);
        downsample@SecmCoords(obj,downrate); 
        obj.image = obj.Ceff.*(obj.image(idxlist,idxlist));
    end

    function lines = line_project(obj,params)
    % lines = obj.LINE_PROJECT(params) Line scan of SECM image.
    % lines = obj.LINE_PROJECT() Line scan of with clpconfig parameters. 
    %
    % The algorithm consists of four steps:
    % 1. Scan with CLP of slope (cos(t),sin(t)) of input angle t.
    %    After rotate the image by input angles, each scan starts from the
    %    largest y-index to the smallest.
    % 2. Shift each scanlines with input shifts.
    % 3. Entrywise multiply each scanlines with input intensity.
    % 4. Convolve each scanlines with input point spread function.
    %
    % LINE_PROJECT has an adjoint operator called back projection.
    %
    % See also SCANLINES.BACK_PROJECT, CLPCONFIG
        if nargin == 1
            params = ProbeParams('config',obj.ticks);
        end
        angles = params.angles.value;
        nlines = length(angles);
        currents = zeros(obj.nmeasures, nlines);

        % Line projection with assigned slope angles
        for I = 1:nlines
            currents(:,I) = flip(sum(obj.rotate_image(obj.image,-angles(I)),2));                     
        end
        lines = ScanLines(obj.ticks, currents, params);

        % Shift each lines with assigned shifts
        if params.shifts.isactive
            s = params.shifts.func(params.shifts.value);
            lines.shift(s);
        end

        % Multiply the current data with specified intensity
        if params.intensity.isactive
            m = params.intensity.func(params.intensity.value);
            lines.mult(m);
        end

        % Convolve the point spread function
        if params.psf.isactive
            p = params.psf.func(params.psf.value);
            lines.conv(p);
            if norm(p) == 0
                error('Point spread function should be nonzero');
            end
        end  
    end

    function draw_image(obj)
    % obj.DRAW_IMAGE(); Draws image with coordinate system.
        h = gca;
        x = obj.ticks;
        y = obj.ticks;
        imagesc(x,y,obj.image);
        set(h,'YDir','normal');
        xlabel('Distance/mm'); xtickformat('%.1f');
        ylabel('Distance/mm'); ytickformat('%.1f');
        c = colorbar(); c.Label.String = 'Currents/A';
                        c.FontSize = h.FontSize - 4;
                        c.Label.FontSize = h.FontSize;
    end
end

methods (Static)
    function img = rotate_image(img,angle)
    % img = obj.ROTATE_IMAGE(img,angle) Rotate image counterclockwise.
        img = fft_rotate(img,angle);
    end
end

methods (Access = private)
    function check_issame(obj1,obj2)
        if ~isequal(obj1.ticks,obj2.ticks)
            error('SecmImage convolution available only in same coordinates')
        end
    end
end

end % classdef