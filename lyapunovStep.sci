// Lyapunov Spectrum
// Neighborhood: neighborhood e_Radius/ , ref ,as many neighbors as dimension.
// one radius per neighborhood (step)

function [jacobian] = lyapunovStep(ps,neighborhood,lag)
    l = length(neighborhood.NNeighbors);
    V = (ps(:,neighborhood.NNeighbors) - ps(:,neighborhood.Ref).*.ones(1,l))';
    s = size(ps);
    if rank(V) < s(1) then
        disp(rank(V))
    end
    lin_map = (ps($,neighborhood.NNeighbors + lag*ones(1,l)) - ps($,neighborhood.Ref+lag).*.ones(1,l))';
    a = inv(V'*V)*V'*lin_map;

    jacobian = [zeros(neighborhood.dim-1,1),eye(neighborhood.dim-1,neighborhood.dim-1);a'];
    
endfunction
