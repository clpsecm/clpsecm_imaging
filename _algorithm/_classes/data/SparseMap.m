classdef SparseMap < SecmImage
% SPARSEMAP Creates sparse map under clp-secm imaging system.
%   obj = SPARSEMAP(ticks,'random-location',radius,ndicts) Creates the map uniform
%   randomly with given dictionary radius and numbers.
% 
%   obj = SPARSEMAP(ticks,'random-all',radius,ndicts) Creates the map uniform
%   randomly with given dictionary radius and numbers and 
%   magnitude~N(1,1/2).
% 
%   obj = SPARSEMAP(ticks,'userdef',locations,magnitude) Create map by its 
%   input location and magnitude infos.
%
% SPARSEMAP methods:
%   SPY_MAP   - Show sparse pattern of map
%   
% SPARSEMAP public fields:
%   MAGNITUDE - magnitude of each dicts
%   LOCATIONS - (x,y) of dict locations
%   NDICTS    - number of dict
%
% For other methods, see SecmImage
%   
% See also SECMIMAGE
properties
    magnitude % vector(ndict,1); magnitude of each dicts
    locations % vector(ndict,2); location of each dicts
    ndicts  % scalar; number of dicts
end

methods
    function obj = SparseMap(ticks, varargin)
        % SPARSEMAP Construct the sparsity map
        gentype = varargin{1};
        if ~strcmp(gentype,'random-location') && ...
           ~strcmp(gentype,'random-all') && ...
           ~strcmp(gentype,'userdef')
            error('Wrong map generation type.');
        end

        % Assign input parameters
        obj   = obj@SecmImage(ticks);
        image = zeros(obj.nmeasures, obj.nmeasures);
        xylim = max(ticks);
        if strcmp(gentype,'random-location')
            radius = varargin{2};
            ndicts = varargin{3}; 
            magnitude = ones(ndicts,2);
        elseif strcmp(gentype,'random-all')
            radius = varargin{2};
            ndicts = varargin{3};
            magnitude = 1+randn(ndicts,1)/4;
        elseif strcmp(gentype,'userdef')
            locations = varargin{2};
            magnitude = varargin{3};
            ndicts = size(locations,1);
        end

        % Unif-random / Self assign the sparsity map
        idict = 1;
        while idict <= ndicts
            % Find location
            if strcmp(gentype,'random-location') || ...
               strcmp(gentype,'random-all') 
                r = (xylim-2*radius)*sqrt(rand());
                t = 2*pi*rand();
                [xloc,xidx] = closest(ticks,r*cos(t));
                [yloc,yidx] = closest(ticks,r*sin(t));
            elseif strcmp(gentype,'userdef')
                [xloc,xidx] = closest(ticks,locations(idict,1));
                [yloc,yidx] = closest(ticks,locations(idict,2));
            end
            
            % Update map
            if image(yidx,xidx) == 0  % The location not assigned before
                image(yidx,xidx) = magnitude(idict);
                locations(idict,:) = [xloc,yloc];
                idict = idict+1;
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

