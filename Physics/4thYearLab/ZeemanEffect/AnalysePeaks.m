function [averageDeltaV, averageErrorDeltaV] = AnalysePeaks(peakCenter, peakErrors, filename)

fabryPerotD = 3/1000;
refIndex = 1.456;
errorSmallDelta = double.empty(4, 0);

%Calculating errors from center of peaks
firstOrderRadii = [abs((peakCenter(4,1) - peakCenter(5,1))/2), abs((peakCenter(4,2) - peakCenter(5,2))/2)];
secondOrderRadii = [abs((peakCenter(3,1) - peakCenter(6,1))/2), abs((peakCenter(3,2) - peakCenter(6,2))/2)];
thirdOrderRadii = [abs((peakCenter(2,1) - peakCenter(7,1))/2), abs((peakCenter(2,2) - peakCenter(7,2))/2)];
fourthOrderRadii = [abs((peakCenter(1,1) - peakCenter(8,1))/2), abs((peakCenter(1,2) - peakCenter(8,2))/2)];

% firstOrderRadii = [abs((1024 - peakCenter(4,1))), abs((1024 - peakCenter(4,2)))];
% secondOrderRadii = [abs((1024 -peakCenter(3,1))), abs((1024 - peakCenter(3,2)))];
% thirdOrderRadii = [abs((1024 - peakCenter(2,1))) , abs((1024 - peakCenter(2,2)))];
% fourthOrderRadii = [abs((1024 - peakCenter(1,1))) , abs((1024 - peakCenter(1,2)))];
% 


%Errors on radius, by propagating erroerror1rs from peak centers
firstOrderRadiiError = [sqrt((peakErrors(4,1)^2 + peakErrors(5,1)^2 )), sqrt((peakErrors(4,2)^2 + peakErrors(5,2)^2 ))];
secondOrderRadiiError = [sqrt((peakErrors(3,1)^2 + peakErrors(6,1)^2 )), sqrt((peakErrors(3,2)^2 + peakErrors(6,2)^2 )) ];
thirdOrderRadiiError = [sqrt((peakErrors(2,1)^2 + peakErrors(7,1)^2 )), sqrt((peakErrors(2,2)^2 + peakErrors(7,2)^2 ))];
fourthOrderRadiiError = [sqrt((peakErrors(1,1)^2 + peakErrors(8,1)^2 )), sqrt((peakErrors(1,2)^2 + peakErrors(8,2)^2 ))];


%calculating error on deltas
radii = [firstOrderRadii', secondOrderRadii', thirdOrderRadii', fourthOrderRadii' ];
errors = [firstOrderRadiiError', secondOrderRadiiError', thirdOrderRadiiError', fourthOrderRadiiError'];





delta_1a = secondOrderRadii(1)^2 - firstOrderRadii(1)^2;
delta_1b = secondOrderRadii(2)^2 - firstOrderRadii(2)^2;

delta_2a = fourthOrderRadii(1)^2 - thirdOrderRadii(1)^2;
delta_2b = fourthOrderRadii(2)^2 - thirdOrderRadii(2)^2;

delta_1 = mean([delta_1a, delta_1b]);
delta_2 = mean([delta_2a, delta_2b]);

averageDelta = mean([delta_1, delta_2]);
[errorSmallDelta, errorDelta] = deltaErrors(radii, errors);

DeltaK1 = abs(firstOrderRadii(2)^2 - firstOrderRadii(1)^2)/(averageDelta*2*fabryPerotD *refIndex);
DeltaK2 = abs(secondOrderRadii(2)^2 - secondOrderRadii(1)^2)/(averageDelta*2*fabryPerotD *refIndex);
DeltaK3 = abs(thirdOrderRadii(2)^2 - thirdOrderRadii(1)^2)/(averageDelta*2*fabryPerotD *refIndex);
DeltaK4 = abs(fourthOrderRadii(2)^2 - fourthOrderRadii(1)^2)/(averageDelta*2*fabryPerotD *refIndex);


errorDeltaK1 = sqrt(((errorSmallDelta(1))/(2 * refIndex * fabryPerotD * averageDelta))^2 + ((errorDelta * abs(firstOrderRadii(2)^2 - firstOrderRadii(1)^2) )/(2 * refIndex * fabryPerotD * averageDelta^2))^2);
errorDeltaK2 = sqrt(((errorSmallDelta(2))/(2 * refIndex * fabryPerotD * averageDelta))^2 + ((errorDelta * abs(secondOrderRadii(2)^2 - secondOrderRadii(1)^2) )/(2 * refIndex * fabryPerotD * averageDelta^2))^2);
errorDeltaK3 = sqrt(((errorSmallDelta(3))/(2 * refIndex * fabryPerotD * averageDelta))^2 + ((errorDelta * abs(thirdOrderRadii(2)^2 - thirdOrderRadii(1)^2) )/(2 * refIndex * fabryPerotD * averageDelta^2))^2);
errorDeltaK4 = sqrt(((errorSmallDelta(4))/(2 * refIndex * fabryPerotD * averageDelta))^2 + ((errorDelta * abs(fourthOrderRadii(2)^2 - fourthOrderRadii(1)^2) )/(2 * refIndex * fabryPerotD * averageDelta^2))^2);

averageDeltaK = mean([DeltaK1, DeltaK2, DeltaK3, DeltaK4]);
averageErrorDeltaK = sqrt(sum([errorDeltaK1^2, errorDeltaK2^2, errorDeltaK3^2, errorDeltaK4^2]))/4;
averageDeltaV = averageDeltaK %* (3*10^8);
averageErrorDeltaV = averageErrorDeltaK %* (3*10^8);

end

