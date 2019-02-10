function [values, sliceDictionary, N] = gaussFit(x, y, peaks)

fringeDictionary = containers.Map();
sliceDictionary = containers.Map();

% gaussFit = double.empty();

ringIndex = 1;
peakIndex = 2 ;
sliceIndex = 1;
N=double.empty();
while 1
    try

        ringEnds(ringIndex) = (peaks(peakIndex) + peaks(peakIndex+1)) / 2;

        
        xValues = x(sliceIndex:round(ringEnds(ringIndex)));
        yValues = y(sliceIndex:round(ringEnds(ringIndex)));
%         N =  pointSpread(xValues, yValues)
        N = [N, pointSpread(xValues, yValues)];
        sliceIndex = round(ringEnds(ringIndex));
        
%         if ringIndex == 1
%             newX = xValues;
%             newY = yValues;
%         end
        
    
        f = fit(xValues, yValues, 'gauss2');
        fringeDictionary(strcat('peak',num2str(ringIndex))) = f;
        sliceDictionary(strcat('peak',num2str(ringIndex))) = xValues;

        
        ringIndex = ringIndex + 1;
        peakIndex = peakIndex + 2;
        
%         fprintf(num2str(ringIndex));
        
    catch
        break
    end 
    
end
values = fringeDictionary;

end