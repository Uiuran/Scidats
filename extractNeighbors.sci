function  [neighbors] = extractNeighbors(phspace,nbhood)
    
    l = length(nbhood.NNeighbors);
    
    neighbors = zeros(nbhood.dim,l);
    
    for i = 1:l
       
       neighbors(:,i) = phspace(1:nbhood.dim,nbhood.NNeighbors(i));
       
    end    
    
endfunction
