//    Update a neighborhood by changing its e_radius or its dimensions and update the list of NNeighbors
function [neighborhood0]= updateNeighborhood(ps,neighborhood,field,value,varargin)
    
    [lhs,rhs] = argn();
    
    neighborhood0 = neighborhood;
    pss = size(ps);
    select field
    case "e_Radius" then
        neighborhood0.e_Radius = value;
        neighborhood0.NNeighbors = find(neighborhood.Dlist < value & neighborhood0.Dlist ~= 0 );
                    
        if rhs == 5 then    
        iso = isochron(neighborhood0.Ref,neighborhood0.NNeighbors,varargin(1));
        else
        iso = isochron(neighborhood0.Ref,neighborhood0.NNeighbors,10);
    end
    
    
    if rhs == 6 then
        lag = varargin(2);
        l = length(neighborhood0.NNeighbors);
        iso = [iso,find(pss(2)*ones(1,l) - neighborhood0.NNeighbors - lag*ones(1,l) < 0),];
    end
    
        neighborhood0.NNeighbors(iso) = [];    
            
    case "dim" then
        neighborhood0.dim = value;
        if neighborhood0.dim == neighborhood.dim + 1  then
        neighborhood0.Dlist = sqrt(dlist(ps(1:value,:),'euclid',[1,length(ps(1,:))],neighborhood0.Ref,value,1,neighborhood.Dlist.*neighborhood.Dlist));
        neighborhood0.NNeighbors = find(neighborhood0.Dlist < neighborhood0.e_Radius & neighborhood0.Dlist ~= 0);
       
        if rhs == 5 then    
        iso = isochron(neighborhood0.Ref,neighborhood0.NNeighbors,varargin(1));
        else
        iso = isochron(neighborhood0.Ref,neighborhood0.NNeighbors,10);
    end
    
        
    if rhs == 6 then
        lag = varargin(2);
        l = length(neighborhood0.NNeighbors);
        iso = [iso,find(pss(2)*ones(1,l) - neighborhood0.NNeighbors - lag*ones(1,l) < 0),];
    end
    
    
        neighborhood0.NNeighbors(iso) = [];
    else
        error('The dimension of the updated neighborhood must be incremented by one in relation to the previous.');
    end
    
    else   
    error("Field not found");
        
        
    end
    
endfunction
