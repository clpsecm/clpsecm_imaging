classdef DictProfile < SecmImage
% DICTPROFILE Creates dictionary profile under clp-secm imaging system.
%   obj = DICTPROFILE(ticks,func) Creates a dictionary element using input
%   funtion handle func under coordinate system depend on ticks.
%
%   Please check its method in object SecmImage
%
% See also SECMIMAGE   
properties 
    func
end

methods
    function obj = DictProfile(ticks,func)
        % Creat dictionries with provided function handle.
        obj = obj@SecmImage(ticks);
        obj.func = func;
        obj.func_to_image();
    end
    function downsample(obj,downrate)
        % Downsample the image of dictionary 
        downsample@SecmCoords(obj,downrate);
        obj.func_to_image();
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
end

methods (Access = private)
    function func_to_image(obj)
        % Draw the indicator function as a map using monte-carlo method 
        N = obj.nmeasures;
        image = zeros(N,N);
        ns = 50;
        for ix = 1:N
            for iy = 1:N
                x = obj.ticks(ix); 
                y = obj.ticks(iy);
                xs = obj.resolution * (rand(ns,1)-0.5) + x;
                ys = obj.resolution * (rand(ns,1)-0.5) + y;
                image(iy,ix) = sum(obj.func(xs,ys))/ns;
            end
        end
        obj.image = image;
    end
end
end