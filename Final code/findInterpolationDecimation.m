%{
findInterpolationDecimation.m
Autor: Laurens Le Jeune and Jonathan Luijsmans
%}

function [I,D] = findInterpolationDecimation(proportion, accuracy)

    proportion = round(proportion,accuracy);
    temporaryI = 1;
    while mod(proportion,floor(proportion)) ~= 0 && ~(temporaryI == 10^accuracy)
        %While there numbers other than zeros behind the comma,
        %keep on multiplying with 10
        %Since doubles are being used, the second argument is necesaary to
        %prevent infinite loops. (Example, a 94 was stored as
        %93.99999999999... which causes and infinite loop
        temporaryI = temporaryI*10;
        proportion = proportion * 10;
    end
    %Find the greatest common divisor:
    GCD = gcd(round(proportion),temporaryI);
    %The values for I and D can then be found
    D = temporaryI / GCD;
    I = int32(proportion / GCD);
end