//falseNearestNeighbor(signal,param_r,test_dims)
//signal: the scalar time series 
//params: the parameter that determines if a certain neighbor is false, a neighbor is false if distance(m+1)/distance(m) > params(2), m is the tested dimension, the params(1) is the radius of the neighborhood where you find the nearest neighbor for the reference, providing this parameter makes the algorithm more efficient.
//test_dims: a integer, last dimension to be tested, if test_dims = 16, the algorithm will test for false neighbors until dimension 16.
//iso_time: a integer, the mean period that a orbit takes to return to a neighborhood, use it to not include isochronous data set in the calculation.
//varargin(1): delay of the tested dimension, optional, if not given will be considered delay = 1.
//varargin(2): a integer nnum that tells the algorithm the aproximated number of neighborhoods to be used (actually the exact number is floor(N/njump) where njump = floor(N/varargin(2))
function [fnn] = falseNearestNeighbor(signal,params,test_dims,iso_time,varargin)
    
    [lhs,rhs] = argn();
    if rhs == 5 then       
        delay = varargin(1);               
        data_std = st_deviation(signal);
        vecset = embedding(signal,test_dims,delay);    
        l = size(vecset);
        N = l(2);
        nnum = N;
    elseif rhs == 6 then        
        delay = varargin(1);
        nnum = varargin(2);
        data_std = st_deviation(signal);
        vecset = embedding(signal,test_dims,delay);    
        l = size(vecset);
        N = l(2);    
    else
        
        delay = varargin(1);               
        data_std = st_deviation(signal);
        vecset = embedding(signal,test_dims,delay);    
        l = size(vecset);
        N = l(2);
        nnum = N; 
        
    end         
    counter = zeros(1,test_dims);    

//    Select a number of nnum equally spaced neighborhoods jumping njump elements
//between the central (reference) points of  the neighborhood
    njump = floor(N/nnum);
    refs = 1:njump:N;
//    The number of references equally spaced by njump may be greater than the given
//nnum, so we must update nnum to the value.
    nnum = length(refs);    

//    Loop through the neighborhoods
//    subset = [];
    [sorted,dim1orig] = gsort(vecset(1,:),'c','i');
    for j = 1:nnum
        uppernear = 0;  
        nearest = 0;    
        // First dimension step
        refsort = find(refs(j) == dim1orig);
        next = refsort + 1;

        near = 1;        
        if next <= N then
                neighbordist = abs(sorted(refsort)-sorted(next));
                while neighbordist < params(1)

                        // Exclude isochronous points
                        if abs(refs(j) - dim1orig(next)) > iso_time then
                           //subset = [subset;next];   obsolet for a while, TODO - use this code to implement saturation of points in neighborhoods, as a method to detect embedding dimension

                            // The conditional on params(1) warrants for too close (or zero distance) points to not compromise the algorithm.,
                            //it also detect the nearest point to the reference
                            if near == 1 & neighbordist > params(1)/1000 then
                                uppernear = next;
                                near = 0;
                                break; // since subset is obselot
                            end    
                        end
                        next = next + 1;
                        if next > N then
                            break;
                        end
                        neighbordist = abs(sorted(refsort)-sorted(next));
                end  
        end    
        
        // Similar search is performed for points in the lower part, relative to reference, of sorted array,
        // it also compare the nearest point of the lower part with the former of the upper part.
        near = 1;
        next = refsort - 1;
        if next >= 1 then
            neighbordist = abs(sorted(refsort)-sorted(next));
            while neighbordist < params(1)

            // Exclude isochronous points
                if abs(refs(j) - dim1orig(next)) > iso_time then                
//                    subset = [next;subset];
                    if near == 1 then                                      
                        if uppernear ~= 0 then                 
                            if (neighbordist <= abs(sorted(refsort)-sorted(uppernear))) & (neighbordist > params(1)/1000) then
                                nearest = next;
                                near = 0;
                                else 
                                nearest = uppernear;
                                near = 0;
                            end
                        elseif neighbordist > params(1)/1000
                            nearest = next;
                            near = 0;
                        end
                    end
                end
                next = next - 1;
                if next == 0 then
                    break;
                end
                neighbordist = abs(sorted(refsort)-sorted(next));
            end
        end
        if nearest == 0 then
            disp(nearest,'Nearest point');            
            disp('No nearest neighbor found in the region of sphere with radius between params(1)/1000 and params(1), you must choose another params(1) and restart, to go ahead');
            fnn = counter/j;
            [fnn,ref] = return(fnn,j);
        end
        // One must copy the neighborhood vectors in order to use subset for dimension 2 and compare with the former.    
//        copyset = subset;
//        subset = [];
        // Further dimensional steps for the reference refs(j)        
        distlist = sqrt( dlist(vecset,'euclid',[1,N],refs(j),2) );          
        nearest = dim1orig( nearest );
        dnn_dimbefore = sqrt(distance(vecset,refs(j),nearest,'euclid',1));
        for i= 2:test_dims        
            
//            l2 = length(subset);            
//            near = 1;
//            for neighb = 1:l2
//                distlist(neighb) = sqrt(distance(vecset,refs(j),dim1orig( subset(neighb) ),'euclid',i)); 
//            if distlist(neighb) < params(1) then
//              subset = [subset;copyset(neighb)];            
//            end
//            end
            // the counter receives a unit if the element of previous dimension is not found in the neighborhood
//            if length(subset) => 1 then a
//            el = find(nearest == subset);
//            else 
//            el = [];
//            endneighbordist = abs(sorted(refsort)-sorted(next));
            next = 1;
//            test = find(subset == nearest);
            if distlist( nearest )/dnn_dimbefore > params(2) then                
                        counter(i-1) = counter(i-1)+1;
            end
       
            // Sort the distance array to find the nearest neighbor in this dimension
            [nearer,norig] = gsort(distlist,'c','i');                     

            // The loop ahead is a warranty that nearest neighbor is not at distance < params(1) from reference            
            nearest = norig(next);
            dnn_dimbefore = nearer(next);
            while dnn_dimbefore < params(1)/1000
                next = next + 1;
                if next > N then
                    break;
                end
                nearest = norig(next);
                dnn_dimbefore = nearer(next);
            end
            if dnn_dimbefore < params(1)/1000 then
                disp('Delay vectors are concentred too near from reference, stopping the calculation');
                fnn = counter/j;
               [fnn,dim,dl,ref] = return(fnn,i,distlist,j);
            end
            
            if i ~= test_dims then
            distlist = sqrt( dlist(vecset,'euclid',[1,N],refs(j),i+1,1,distlist.*distlist) );
            end
            
        end

//        subset = [];
    end
    fnn=counter/nnum;    
    
endfunction
