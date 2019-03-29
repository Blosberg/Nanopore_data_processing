function result = MF(Pr, CF)
    result = (Pr*CF)./(1 + (CF*(Pr-1)) );
end