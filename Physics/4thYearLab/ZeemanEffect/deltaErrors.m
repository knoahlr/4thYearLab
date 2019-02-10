function [errorSmallDelta, errorDelta] = deltaErrors(radii, errors)

errorSmallDelta1 = sqrt( (2*radii(2,1)*errors(2,1))^2 +  (2*radii(1,1)*errors(1,1))^2 );
errorSmallDelta2 = sqrt( (2*radii(2,2)*errors(2,2))^2 +  (2*radii(1,2)*errors(1,2))^2);
errorSmallDelta3 = sqrt( (2*radii(2,3)*errors(2,3))^2 +  (2*radii(1,3)*errors(1,3))^2);
errorSmallDelta4 = sqrt( (2*radii(2,4)*errors(2,4))^2 +  (2*radii(1,4)*errors(1,4))^2);

errorSmallDelta = [errorSmallDelta1, errorSmallDelta2, errorSmallDelta3, errorSmallDelta4];


errorDeltaB1 = sqrt( (2*radii(2,2)*errors(2,2))^2 +  (2*radii(2,1)*errors(2,1))^2 );
errorDeltaB2 = sqrt( (2*radii(2,4)*errors(2,4))^2 +  (2*radii(2,3)*errors(2,3))^2 );

errorDeltaA1 = sqrt( (2*radii(1,2)*errors(1,2))^2 +  (2*radii(1,1)*errors(1,1))^2 );
errorDeltaA2 = sqrt( (2*radii(1,4)*errors(1,4))^2 +  (2*radii(1,3)*errors(1,3))^2 );

errorDeltaB = sqrt(sum([errorDeltaB1^2, errorDeltaB2^2]))/2;
errorDeltaA = sqrt(sum([errorDeltaA1^2, errorDeltaA2^2]))/2;

errorDelta = sqrt(sum([errorDeltaB^2, errorDeltaA^2]))/2;


end

