//     Calculate the mean time interval between maxima and minima  of the signal by making use of its derivative
// and histograming around the time intervals. The user must provide 3 different bin numbers, to the algorithm
// calculate the mean value of the mean time interval, of three estimated histograms

function [meantime,time] = meanPeriod(signal,signalderivative,binnum,distminimal,treshold)
   
    if length(binnum) ~= 3 then
        error('wrong number of bin numbers for histograming, 3 expected');        
    end
    meantime = 0;
[orderedsd,time] = gsort(abs(signalderivative));
time = gsort(time($-treshold:$));

index = 1;
l = length(time);

while index < l

d = time(index)-time(index+1);

if d < distminimal then
time(index+1) = [];
l = length(time);

else
index = index + 1;

end
end

d = diff(time);

d = gsort(abs(d));

hist0 = [0;0;0];
hist1 = [0;0;0];
hist2 = [0;0;0];
binsize0 = floor((d(1)-d($))/binnum(1));
binsize1 = floor((d(1)-d($))/binnum(2));
binsize2 = floor((d(1)-d($))/binnum(3));

l = length(d);
index = 1;
axis0 = d($);
axis1 = d($);
i = 0;
while i ~= l
    
    axis1 = axis1 + binsize0;
    hist0(1,index) = 0;
    hist0(2,index) = axis0;
    hist0(3,index) = axis1;    
    
    while d($-i) < axis1 & axis0 <= d($-i)
        i = i + 1;
        hist0(1,index)= hist0(1,index) + 1;
        if i == l then
            break;
        end
    end
    hist0(1,index) = hist0(1,index)/l;
    index = index + 1;
    axis0 = axis0 + binsize0;
                
end

   
index = 1;
axis0 = d($);
axis1 = d($);
i = 0;
while i ~= l
    
    axis1 = axis1 + binsize1;
    hist1(1,index) = 0;
    hist1(2,index) = axis0;
    hist1(3,index) = axis1;
    
    while d($-i) < axis1 & axis0 <= d($-i)
        i = i + 1;
        hist1(1,index)= hist1(1,index) + 1;
        
        if i == l then
            break;
        end
    end
    hist1(1,index) = hist1(1,index)/l;
    index = index + 1;
    axis0 = axis0 + binsize1;
                
end

index = 1;
axis0 = d($);
axis1 = d($);
i = 0;
while i ~= l
    
    axis1 = axis1 + binsize2;
    hist2(1,index) = 0;
    hist2(2,index) = axis0;
    hist2(3,index) = axis1;
    while d($-i) < axis1 & axis0 <= d($-i)
        i = i + 1;
        hist2(1,index)= hist2(1,index) + 1;
        if i == l then
            break;
        end
    end
    hist2(1,index) = hist2(1,index)/l;
    index = index + 1;
    axis0 = axis0 + binsize2;
                
end

[hord0] = gsort(hist0,'lc');
[hord1] = gsort(hist1,'lc');
[hord2] = gsort(hist2,'lc');

m0 = hord0(2,1) + floor((hord0(3,1) - hord0(2,1))/2);
m1 = hord1(2,1) + floor((hord1(3,1) - hord1(2,1))/2);
m2 = hord2(2,1) + floor((hord2(3,1) - hord2(2,1))/2);

meantime = floor((m0+m1+m2)/3);
   
endfunction
