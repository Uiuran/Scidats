function [cmass] = centerofmass(nneighbors)
   
   l = size(nneighbors);
   cmass = zeros(1:l(1))';

   for w=1:l(2)

       for j=1:l(1)
       cmass(j) = cmass(j) + nneighbors(j,w);
       end
       
   end

   cmass = cmass/l(2);

endfunction
