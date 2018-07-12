function newvar = deepcopy(var)
% newvar = DEEPCOPY(var) Copy variable when it's may or may not copyable.
if iscell(var)
    newvar = cell(size(var));
    for I = 1:numel(var)
        try newvar{I} = copy(var{I}); 
        catch; newvar{I} = var{I}; end
    end
else
    try newvar = copy(var); 
    catch; newvar = var; end
end


