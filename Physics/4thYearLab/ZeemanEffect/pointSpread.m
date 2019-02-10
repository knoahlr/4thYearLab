function N = pointSpread(x,y)


threshhold = 10;


for j=length(y):-1:1
    if y(j) > threshhold 
        xRightLim = x(j);
        break;
    end
end
for i=1:length(y)
    if y(i) > threshhold
        xLeftLim = x(i);   
        break;
    end
end


N = xRightLim - xLeftLim;
end
