function [vec] = embedding(time_series,veclength,lag)
    i = 0;
    k = 0;
    vcont = 0;
    lenght = size(time_series);
        
        while i < lenght(1) - lag*(veclength-1)
            vcont = vcont + 1;
            k = i;
            
                while ((k-i+1) < (veclength + 1)) & ((i+(k-i)*lag) < lenght(1))
                    vec(k-i+1,vcont)= time_series(i+1+(k-i)*lag);
                    k = k + 1;
                end
            i = i + 1;
        end
endfunction
