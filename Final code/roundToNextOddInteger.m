%{
roundToNextOddInteger.m
Autor: Laurens Le Jeune and Jonathan Luijsmans
%}

function m = roundToNextOddInteger(n)
%This function rounds up a given number n. 
    n = ceil(n);
%If the result is even, it needs to be incremented.
    if(mod(n,2) == 0)
        m = n + 1;
    else
        m = n;
    end
end