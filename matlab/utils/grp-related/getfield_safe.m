function val = getfield_safe(opt,field,default_val,debug_str)
%
% Get the field item if exists, and return default otherwise.
%
global VERBOSE
if ~exist('VERBOSE','var')
    VERBOSE = 0;
end

if isfield(opt,field) % if field exists, use it
    val = getfield(opt,field);
else % otherwise, use the default value
    if VERBOSE
        fprintf(2,'[%s][%s] field does not exist. Using [%s] instead. \n',...
            debug_str,field,num2str(default_val));
    end
    val = default_val;
end
