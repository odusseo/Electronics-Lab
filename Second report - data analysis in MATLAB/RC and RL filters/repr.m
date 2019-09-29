% str = repr(var)
%
% Converts var to a string representation str.
%
% This code is inherited from UTN, mainly from D Nicolodi
%
function str = repr(var)
  
  if ischar(var)
    str = var;
    return
  end
  
  if isscalar(var)
    str = sprintf('%g',var);
    return
  end
  
  if isnumeric(var)
    str = mat2str(var);
    return
  end
  
  str = '';
  
end