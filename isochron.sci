function [iso] = isochron(ref,neighborhood,isotime)

 iso = find(abs(ref - neighborhood) < isotime);

endfunction
