function nh = HarmIndexVector(nH) 

counter = 1;
nh(counter) = 0;
counter = counter+1;
coeff = 1;
for i = 1:nH
    nh(counter) = coeff;
    nh(counter+1) = coeff;
    coeff = coeff+1;
    counter = counter+2;
end

end