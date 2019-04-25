classdef IPalmSecmRealdata < IPalm
    properties
        D  % true kernel
        lines % observed lines
    end
    
    methods
        function obj = IPalmSecmRealdata(problem,disc,lines)
            obj = obj@IPalm(problem);
            obj.D = disc;
            obj.lines = lines;
        end

        function display_result(obj)
            subplot(221);
            display_result@IPalm(obj);
            title('Objectivel value');
            
            subplot(222);
            Yhat = obj.vars{1}*obj.D;
            Yhat.draw_image();
            title('Current image Y');
            
            subplot(223);
            obj.lines.plot_psf();
            title('Point spread function');
            legend off;
            
            subplot(224);
            obj.vars{1}.draw_image; 
            set(gca,'YDir','normal');
            title('Current map X0');
            
            angles = obj.vars{2}.angles.value;
            if iscolumn(angles); angles = angles'; end
            disp(['angles(deg): ', num2str(angles,'%.2f  ')]);
            
            shifts = obj.vars{2}.shifts.value;
            if iscolumn(shifts); shifts = shifts'; end
            disp(['shifts(mm) : ', num2str(shifts,'%.2f  ')]);
            
            intensity = obj.vars{2}.intensity.value;
            if iscolumn(intensity); intensity = intensity'; end
            disp(['intensity  : ', num2str(intensity,'%.2f  ') ]);
            
            psf = obj.vars{2}.psf.value;
            if iscolumn(psf); psf = psf'; end
            disp(['psf        : ', num2str(psf,'%.2f  ') ]);
        end
    end
end

