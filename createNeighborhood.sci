//    Create a neighborhood, in relation to a reference vector ref pertained to phase space pspace, with dimension 1.
//Note that a neighborhood is always created with one dimension (first component of each vector), to expand the neighborhood
//to high-dimensions one must do updateNeighborhood.The algorithm attemps to exclude the isochronous(to the reference) vectors
//, i.e. vectors that are nearly tangencial to the reference orbit, by projecting the relative(to the reference) vectors into 
// the first forward and backward isochronous relative vectors. The isochronous exclusion is made in order to avoid acounting 
//for properties that depends only of a orbit, and is not a convergent property of the attractor, e.g. Lyap Exponent of the tangencial
// space is always 0.
//    varargin(1) = integer number, if defined, exclude vectors that are less than varargin(1) time units near from reference, can be used with
// meanPeriod function to exclude isochronous vectors in order to evaluate orbit independent property of the attractor (e.g. lyapunov exponents).
//    varargin(2) = type of the distance, calculated for 1 dimension in the list dlist, of the various vectors relative to the reference, can be 
// 'euclid' for euclidian or 'max' for supremum or infinity distance.
// pspace = struct("Vector_set",[],"m",16,"time_delay",3)
//    varargin(3) = delay choosen to build coordinated vectors, since you need that the neighbors have at least one plus lag isochrounous point to do some operations

function [neighborhood] = createNeighborhood(pspace,ref,e_rad,varargin)
  
    [lhs,rhs] = argn();
    
    if rhs == 5 then
        typ = varargin(2);
    else
        typ = 'euclid';
    end
    
    pss = size(pspace);
    neighborhood = struct("Neighborhood_index",ref,"Ref",ref,"e_Radius",e_rad,"dim",1,"NNeighbors",[],"Dlist", sqrt(dlist(pspace(1,:),typ,[1,pss(2)],ref,1)) );

    neighborhood.NNeighbors = find(neighborhood.Dlist < neighborhood.e_Radius & neighborhood.Dlist ~= 0.0);    
    
    if rhs >= 4 then
    iso = isochron(ref,neighborhood.NNeighbors,varargin(1));
    else
    iso = isochron(ref,neighborhood.NNeighbors,10);
end

if rhs == 6 then
    lag = varargin(3);
    l = length(neighborhood.NNeighbors);
    iso = [iso,find(pss(2)*ones(1,l) - neighborhood.NNeighbors - lag*ones(1,l) < 0),];
end
    
    neighborhood.NNeighbors(iso) = [];    

endfunction

