//     Calculate the squared distance between two vector, ref and nn, in the dimension dim,
//optionally given as 5th argument. One choose between 'euclid' or 'max' norm in type, By making inner prod of vector differences
//, by recursive calculation in the former case, if the 6th argument is 1, using a pre-calculed
//distance, on a dim minus 1 dimensional space, given as the 7th argument.
//the syntax is:
//
// d = distance(phase_space,reference,neighbortype,type,dimension,recursive,predistance);
//
//     The last tree arguments are optional
function [d]= distance(ps,ref,nn,typ,varargin);
    [um,io] = argn();

    if io == 5 then
        
        dim = floor(varargin(1));
        recursive = 0;
        predist = 0;
        if dim > 0 then        
        
        if typ == 'max' then            
        d = max(abs(ps(1:dim,nn)-ps(1:dim,ref)));
        elseif typ == 'euclid' then
        d =  (ps(1:dim,ref) - ps(1:dim,nn))'*(ps(1:dim,ref) - ps(1:dim,nn));
        end
        else d = (ps(:,ref) - ps(:,nn))'*(ps(:,ref) - ps(:,nn)) ;
        
        end
        
    elseif io < 5 then
        dim = 0;
        recursive  = 0;
        predist = 0;
        if typ == 'max' then            
        d = max(abs(ps(:,nn)-ps(:,ref)));
        elseif typ == 'euclid' then
        d =  (ps(:,ref) - ps(:,nn))'*(ps(:,ref) - ps(:,nn));
    end
    
    else        
        if io == 6 then
            error('The distance, in the dimension-1 space, between two vectors must be given as last argument');
        end
      
        dim = floor(varargin(1));
        recursive = varargin(2);
        predist = varargin(3);
        
        if dim > 0 & recursive == 0 then        
       
        if typ == 'max' then            
        d = max(abs(ps(1:dim,nn)-ps(1:dim,ref)));
        elseif typ == 'euclid' then
        d =  (ps(1:dim,ref) - ps(1:dim,nn))'*(ps(1:dim,ref) - ps(1:dim,nn));
        end
        
        elseif dim == 0 & recursive == 0 then             
        
        if typ == 'max' then
        d = max(abs(ps(:,nn)-ps(:,ref)));
        elseif typ == 'euclid' then
        d =  (ps(:,ref) - ps(:,nn))'*(ps(:,ref) - ps(:,nn));
        end
    
        elseif dim > 0 & recursive == 1 then        
        
        if typ == 'max' then  
        d = max([predist,abs(ps(dim,nn)-ps(dim,ref))]);
        elseif typ == 'euclid' then        
        d =  predist + (ps(dim,ref) - ps(dim,nn))'*(ps(dim,ref) - ps(dim,nn));
        end
    
        elseif dim == 0 & recursive == 1 then 
        
        if typ == 'max' then  
        d = max([predist,abs(ps($,nn)-ps($,ref))]);
        elseif typ == 'euclid' then        
        d =  predist + (ps($,ref) - ps($,nn))'*(ps($,ref) - ps($,nn));
        end
    
        else error('You must set at least 3 arguments');
        end

    end
        
endfunction

// varargin(1) = [index_first_vector,index_last_vector] set vectors inside ps to calculate distances
// varargin(2) =  ref0 Integer between 1 <= ref0 <= ps_size, if defined calculates only the distances relative to reference ref0.
// varargin(3) = dim Integer  between 1<= dim <= ps_size(1), dimension of the vectors whose distance will be calculated
// varargin(4) = 1 or 0, recursive mode, calculates distances by adding a dimension to squared distance value of previous dimension 
// varargin(5) = list of positive reals, the values of squared distances, in a dimension = dim - 1, between the calculated points, must have the same size of the number of distances
// dl = dlist(ps) calculates all possible distances, return a list with upper triangular distance matrix, to recover a element d_ij with j>i, one can use take(dl,i,j,size(dl))
function [dl]= dlist(ps,typ,varargin)
    
    [rhs,lhs] = argn();
    ps_size = size(ps);
    index = 1;
    
    if lhs == 3 then
    vec_interval = varargin(1);
    vec_size = vec_interval(2)-vec_interval(1)+1;
    dl = zeros(1,vec_size*vec_size/2-vec_size);
    
    for ref = vec_interval(1):vec_interval(2)-1       
    for i=ref+1:vec_interval(2)
        
        dl(index) = distance(ps,ref,i,typ);
        index = index + 1;
    end
    end

    elseif lhs == 4 then
    
    vec_interval = varargin(1);
    vec_size = vec_interval(2)-vec_interval(1)+1;    
    ref0 = varargin(2);
    if ref0 > vec_interval(2) | ref0 < vec_interval(1) then
        error('Reference vector outside the interval');
    end
    dl = zeros(1,vec_size);    
    
    for ref = vec_interval(1):vec_interval(2)  
           
        dl(index) = distance(ps,ref0,ref,typ);
        index = index + 1;    
    end
    
elseif lhs == 5 then
    
    dim = varargin(3);
    vec_interval = varargin(1);
    vec_size = vec_interval(2)-vec_interval(1)+1;    
    ref0 = varargin(2);  
    
    if ref0 <= vec_interval(2) & ref0 >= vec_interval(1) then
    dl = zeros(1,vec_size);        
    
    for ref = vec_interval(1):vec_interval(2)  
           
        dl(index) = distance(ps,ref0,ref,typ,dim);
        index = index + 1;    
    end
    
    elseif ref0 == 0 then

    dl = zeros(1,vec_size*vec_size/2-vec_size);
    
    for ref = vec_interval(1):vec_interval(2)-1       
    for i=ref+1:vec_interval(2)
        
        dl(index) = distance(ps,ref,i,typ,dim);
        index = index + 1;
    end
    end

    end
    
elseif lhs == 7 then
    
    dim = varargin(3);
    lista = varargin(5);
    vec_interval = varargin(1);
    vec_size = vec_interval(2)-vec_interval(1)+1;    
    ref0 = varargin(2);  
    
    if ref0 <= vec_interval(2) & ref0 >= vec_interval(1) then
    dl = zeros(1,vec_size);        
    if varargin(4) == 1 then
    
    for ref = vec_interval(1):vec_interval(2)  
           
        dl(index) = distance(ps,ref0,ref,typ,dim,varargin(4),lista(ref-vec_interval(1)+1));
        index = index + 1;    
    end
    elseif varargin(4) == 0 then
    
    for ref = vec_interval(1):vec_interval(2)  
           
        dl(index) = distance(ps,ref0,ref,typ,dim);
        index = index + 1;    
    end    
else error('Recursive mode is 1 or 0');
    
    end
    elseif ref0 == 0 then
    dl = zeros(1,vec_size*vec_size/2-vec_size);
    if varargin(4) == 1 then
        
    for ref = vec_interval(1):vec_interval(2)-1       
    for i=ref+1:vec_interval(2)
        
        dl(index) = distance(ps,ref,i,typ,dim,varargin(4),take(lista,ref,i,vec_size));
        index = index + 1;
    end
end

elseif varargin(4) == 0 then
    
    
    for ref = vec_interval(1):vec_interval(2)-1       
    for i=ref+1:vec_interval(2)
        
        dl(index) = distance(ps,ref,i,typ,dim);
        index = index + 1;
    end
    end

else error('Recursive mode is 1 or 0');    
    
end
end

else
      dl = zeros(1,ps_size(2)*ps_size(2)/2-ps_size(2));
    
    for ref = 1:ps_size(2)-1       
    for i=ref+1:ps_size(2)
        
        dl(index) = distance(ps,ref,i,typ);
        index = index + 1;
    end
    end
    end

endfunction
