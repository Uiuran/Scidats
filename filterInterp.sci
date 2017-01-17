//    Interpolation of a data array by a polynomial, either splines or fourier(fft), and filtering by fft,
// the number of points interpolated and filtered must be both even. If interpolation size (interpsize)
// is of the same size of the original data, than one simple returns the same series filtered(or derivative).
// The user also must provide a windowSize with the size less or equal of the data length.
// One can also specify the infimum and supremum of the data domain in a bidimensional vector xdomain,
// otherwise it will be take as being at [0,1].
//
function [interpdata,firstderiv,xx] = filterInterp(data,windowSize,filtersize,interpsize,interptype,xdomain)   

[lhs,rhs] = argn();
l = length(data);
firstderiv = 0;
if isdef('xdomain') then
    
    if length(xdomain) == 2 then
        
       
       xx = linspace(xdomain(1),xdomain(2),interpsize);
    end
    
else
    
    xx = linspace(0,1,interpsize);
    xdomain(1) = 0;
    xdomain(2) = 1;
end

if  length(windowSize) < 2 then
    
    if  windowSize(1) > l then
        error("Window is too large");
    else
    fdata = data(1:windowSize(1));
    l = length(fdata);
    m = modulo(l,2);

end

else
    
    if  windowSize(2) > l then
        error("Window is too large");
    elseif windowSize(1) > l | windowSize(1) < 1 then
        error("Invalid start point");
    else        
    fdata = data(windowSize(1):windowSize(2));
    l = length(fdata);
    m = modulo(l,2);
    
end

end

if ~isdef('interptype') then
    interptype = 'f';
end

it = interptype;
    if m == 0 then
       
       if filtersize ~= 0 | it == "f" | it == "fourier" | it == "fft" then
       fdata = fft(fdata);
       end
       up = fdata(1:l/2);
       down = fdata(l/2+1:l);
       dup = fdata(1:l/2);
       ddown = fdata(l/2+1:l);
       scale = interpsize/l;
       up($-filtersize/2+1:$) = 0;
       down(2:1+filtersize/2)= 0;
       dup($-filtersize/2+1:$) = [];    
       ddown(2:1+filtersize/2)= [];
       
   if it == "f" | it == "fourier" | it == "fft" then
       fdata = [up;zeros((interpsize-l)/2,1);down(1:1+filtersize/2);zeros((interpsize-l)/2,1);down(2+filtersize/2:$) ];       
       interpdata = real(fft(scale*fdata,1));
       
       if lhs > 1 then
       s = length(interpdata);
       firstderiv = filterInterp(diffDerivative(interpdata,l),s-1,(s-1)/2,s-1,'f');       
       end
   
   elseif it == "s" | it == "splines" then

       if filtersize ~= 0  then
       fdata = real(fft([up;down],1));
       end
     
       x = linspace(xdomain(1),xdomain(2),l)';
       sp = splin(x,fdata);
       [a,b] = interp(xx,x,fdata,sp);
       interpdata = a';
       firstderiv = b';
       
   else 
       
       if filtersize ~= 0 then
       interpdata = real(fft([up;down],1));
       else interpdata = [up;down];
       end

   end       
       
   else
       
       
       if filtersize ~= 0 | it == "f" | it == "fourier" | it == "fft" then
       fdata = fft(fdata);
       end       
       first = fdata(1);
       fdata(1) = [];
       l = l-1;
       
       up = fdata(1:l/2);
       down = fdata(l/2+1:l);
       scale = (interpsize)/l;
       up($-filtersize/2+1:$) = 0;
       down(1:filtersize/2)= 0;

       if it == "f" | it == "fourier" | it == "fft" then
       fdata = [first;up;zeros((interpsize-l)/2,1);down(1:filtersize/2);zeros((interpsize-l)/2,1);down(1+filtersize/2:$) ];
       interpdata = real(fft(scale*fdata,1));
       
       if lhs > 1 then
       s = length(interpdata);
       firstderiv = filterInterp(diffDerivative(interpdata,l),s,s/2,0,'f');       
       end
       
       elseif it == "s" | it == "splines" then
              
       l = l + 1;     
       
       if filtersize ~= 0 then
       fdata = real(fft([first;up;down],1));
       else
       fdata = [first;up;down];
       end
       
       x = linspace(0,1,l)';
       xx = linspace(0,1,interpsize)';
       sp = splin(x,fdata);
       [a,b] = interp(xx,x,fdata,sp);
       interpdata = a';
       firstderiv = b';

       
    else 
       
       if filtersize ~= 0 then
       interpdata = real(fft([first;up;down],1));
       else interpdata = [first;up;down];
       end
    
    end
   
   end
endfunction
