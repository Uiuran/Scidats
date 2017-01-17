// Simple heaviside step function

function [argout] = stepFunction(arg)
    
    if arg >= 0 then
        
        argout = 1;
        
    else
        
        argout = 0;
        
    end
    
endfunction
