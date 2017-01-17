function [firstderiv] = diffDerivative(ts,Nold)
    l = length(ts);
    
    firstderiv = diff(ts)*(l/Nold);
    
endfunction
