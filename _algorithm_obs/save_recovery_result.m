function save_recovery_result(filename, imageNumber)

if nargin == 2
    figure(imageNumber);
end
savefig(filename);
export_fig('transparent','-pdf',filename);
save(filename);
