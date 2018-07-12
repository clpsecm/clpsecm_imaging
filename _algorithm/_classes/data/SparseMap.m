classdef SparseMap < SecmImage
% SPARSEMAP Creates sparse map under clp-secm imaging system.
%   obj = SPARSEMAP(ticks,'random',radius,ndicts) Creates the map uniform
%   randomly with given dictionary radius and numbers.
% 
%   obj = SPARSEMAP(ticks,'userdef',locations) Create map using location. 
%
% SPARSEMAP methods:
%   SPY_MAP   - Show sparse pattern of map
%   
% SPARSEMAP public fields:
%   LOCATIONS - (x,y) locations of sparse map
%   NDICTS    - number of non-zero entries 
%
% For other methods, see SecmImage
%   
% See also SECMIMAGE
properties
    locations
    ndicts
end

methods
    function obj = SparseMap(ticks, varargin)
        % SPARSEMAP Construct the sparsity map
        gentype = varargin{1};
        if ~strcmp(gentype,'random') && ~strcmp(gentype,'userdef')
            error('Wrong map generation type.');
        end

        % Assign input parameters
        obj   = obj@SecmImage(ticks);
        image = zeros(obj.nmeasures, obj.nmeasures);
        xylim = max(ticks);
        switch nargin 
            case 3  % User-defined
                locations = varargin{2};
                ndicts = size(locations,1);
            case 4  % Random
                radius = varargin{2};
                ndicts = varargin{3}; 
                locations = zeros(ndicts,2);
        end

        % Unif-random / Self assign the sparsity map
        iloc = 1;
        while iloc <= ndicts
            % Find location
            if strcmp(gentype,'random')
                r = (xylim-2*radius)*sqrt(rand());
                t = 2*pi*rand();
                [xloc,xidx] = closest(ticks,r*cos(t));
                [yloc,yidx] = closest(ticks,r*sin(t));
            elseif strcmp(gentype,'userdef')
                [xloc,xidx] = closest(ticks,locations(iloc,1));
                [yloc,yidx] = closest(ticks,locations(iloc,2));
            end

            % Update map
            if image(yidx,xidx) ~= 1 % The location not assigned before
                image(yidx,xidx) = 1;
                locations(iloc,:) = [xloc,yloc];
                iloc = iloc+1;
            elseif strcmp(gentype,'userdef') % The location is assigned
                error('Repeated locations');
            end
        end
        obj.image = image;
        obj.locations = locations;
        obj.ndicts = ndicts;
    end

    function downsample(obj, ~)
        % [TODO] Suggest implementing nearset neighbor interpolation.
        error('Sparsity map cannot be downsampled.');
    end

    function secmImage = pos(obj)
        % Output and secmimage object
        secmImage = SecmImage(obj.ticks,obj.image);
        secmImage.pos();
    end

    function secmImage = soft(obj,lda)
        % Output and secmimage object
        secmImage = SecmImage(obj.ticks,obj.image);
        secmImage.soft(lda);
    end

    function spy_map(obj)
        % obj.SPY_MAP Draws sparsity pattern with coordinate system.
        spy(obj.image);
        delta = floor(obj.nmeasures/10);
        idx = 1:delta:obj.nmeasures;
        set(gca,'Xtick',idx,...
                'Xticklabel',obj.ticks(idx),...
                'Ytick',idx,...
                'Yticklabel',obj.ticks(idx),...
                'YDir','normal');
        xlabel('Distance/mm')
        ylabel('Distance/mm')
    end
end
end

