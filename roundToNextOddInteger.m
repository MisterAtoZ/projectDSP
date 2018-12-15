function m = roundToNextOddInteger(n)
    
    n = ceil(n);
    if(mod(n,2) == 0)
        m = n + 1;
    else
        m = n;
    end
end