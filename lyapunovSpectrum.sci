// Lyapunov Spectrum
// Neighborhood: neighborhood e_Radius/ , ref ,as many neighbors as dimension.
// one radius per neighborhood (step)

// lyapunovSpectrum(ps,lag,isotime,lastref,e_Radius,first_Ref)

// ps := phase space, set of vectors of delay coordinates with the same or bigger dimension than embedding dimension of Takens Theorem
// lag := delay, integer number used to build the delay coordinated vectors from the scalar signal.
// isotime := An estimate of the time necessary to discard isochronous vector in relation to the orbit of the reference point of a neighborhood.
// dimension := an integer, reduced dimension in relation to the embedding dimension, is supposed to be the dimension of the estimated jacobian;
// (optional - varargin(1)) e_Radius := Is the first radius of the neighborhoods to be used in the calculation of the lyapunov spectrum, it can be change subsequently to achieve a better precision/number of neighbors
// (optional - varargin(2)) first_Ref := first reference point to be used in order to build a neighborhood and calculate the first local jacobian, if not set reference 1 will be used
// (optional - varargin(3)) last_Ref := last reference point, must be lower than  N-d(2)-lag+1

function [spect,esterr,neighbnum] = lyapunovSpectrum(tseries,lag,isotime,dimension,varargin)

    [rhs,lhs] = argn();

    //d(1) reduced dimension, d(2) embedding dimension, d(2) = (d(1) - 1)*lag + 1
    d = [0,0];
    d(1) = dimension;
    d(2) = (d(1)-1)*lag + 1;   
    spect = zeros(d(1),1);
    ps = embedding(tseries,d(2),1);
    ps0 = ps(1:lag:$,:);
    l = size(ps);
    if lhs == 7 then
        
        e_rad = varargin(1);
        ref = varargin(2);
        lastref = varargin(3);
    
        if e_rad <= 0 then
            error('wrong input for neighborhood radius');        
        end
        
        if ref < 1 | ref > l(2) then
            error('reference is not in the phase space set');
        end       
        
        elseif lhs == 6 then
            e_rad = varargin(1);            
            ref = varargin(2);
            lastref = l(2) - lag;
            if e_rad <= 0 then
                error('wrong input neighborhood radius');
            end 
        elseif lhs == 5 then
            e_rad = varargin(1);
            ref = 1;
        else
            e_rad = 0.1;
            ref = 1;    
    end
           
    last = lastref;
    if lastref + lag > l(2) then
        
        disp('Cannot map a vector that is more distance than the embedding coordinate set.');
        disp('Setting last reference to be V(dataset_size - lag)');
        last = l(2) - lag;
        
    end
    
    if ref > last then
        
        disp('first reference cannot be greater than the last, ressetting' );
        ref = last;
        
    end
    
    //    The first tangent basis is the canonical basis, so the jacobian maps a vector V(i) to a vector V(i+lag) also evolving the tangent canonical basis, we do the QR-Decomposition in such evolved basis getting a Q ortogonal basis to be evolved again and a R upper triangular with positive diagonal elements which are related to the tax of expansion of the basis, they will be used to calculate the lyapunov exponents.
    
    Q = list( eye(d(1),d(1)) );
    R = list(0);
    
    //Since lyapunov spectrum is a natural measure in ergodic systems, then it can be approximated by suficiently large number of observations
    esterr = [];
    neighbnum = [];
    for i = ref:lag:last
        esterr(1:d(1),$+1) = zeros(d(1),1);
        neighbnum($+1) = 0;
        neighborhood = createNeighborhood(ps,ref,e_rad,isotime,'euclid',lag);

    //  Adjust the dimension of the new neighborhood to be equal to embedding dimension l(1)
        while neighborhood.dim < l(1)
            neighborhood = updateNeighborhood(ps,neighborhood,'dim',neighborhood.dim+1,isotime,lag);     
        end
    //  Adjust the radius e_Radius of the neighborhood to have at least as many neighbors as the number of embedding dimensions, this is done so the
    // jacobian can be full rank and have a good statistical acceptance. 
//        while length(neighborhood.NNeighbors) < 2*d(1)
        while length(neighborhood.NNeighbors) < 30
            neighborhood = updateNeighborhood(ps,neighborhood,'e_Radius',neighborhood.e_Radius + 0.01,isotime,lag);
        end
        neighborhood.dim = d(1);        
        jacobian = lyapunovStep(ps0,neighborhood,lag);
        Q_next = jacobian*Q($);
        [Q($+1),R($+1)] = qr(Q_next);
	    predict = jacobian*ps0(:,ref);
	    esterr(1:d(1),$) = abs(predict-ps0(:,ref+lag));
        neighbnum($) = length(neighborhood.NNeighbors);
        // The native qr method, of Scilab, does negative lines in R to order the diagonal in increasing order, also negativing the corresponding Q column, to avert this we do the loop below.
        for j = 1:d(1)            
            if R($)(j,j) < 0 then
                R($)(j,:) = -R($)(j,:);
                Q($)(:,j) = -Q($)(:,j);            
            end
        end
        ref = ref + lag;
    end
    
    // The looping below calculates all exponents
    norma = [];
    numevo = length(Q);
    spect = zeros(d(1),numevo-1);
    
    for i = 2:numevo
        norma(i-1) = i-1;
        tmp = diag(R(i));
        
        spect(:,$:-1:i-1) = spect(:,$:-1:i-1) + log(tmp)*ones(1,numevo-1-i+2);
        
    end
    
    norma = ones(1,numevo-1)./(norma)';
    
    norma = ones(d(1),1)*norma;
    
    spect = spect.*norma;
    
endfunction
