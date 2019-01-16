addpath("C:\Users\Noah Workstation\Desktop\P_PR_1\repo\attentionNN\lidar\atten\sim\noah_sim\");


%Opening and loading fringe data
fringeFolder = "D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Zeeman Effect\Data\Fringes_2";
DataFolder  = "..\Logs\Data\Fringes_3";

foldernames = dir(DataFolder);
gaussCoeff = double.empty();

fileNo = 1;

for i=1:length(foldernames)
    
    
    folder = strcat(DataFolder,filesep,foldernames(i).name);

    
    fp = strcat(folder, filesep, "peaks.npy");
    try
        peaks = readNPY(strcat(folder, filesep, "peaks.npy")); 
        xValues = readNPY(strcat(folder, filesep, "xValues.npy")); 
        yValues = readNPY(strcat(folder, filesep, "yValues.npy")); 
        
        figure(fileNo);

        plot(xValues, yValues, 'b');
        hold on;
        [coeff, slices] = gaussFit(xValues, yValues, peaks);

%         gaussCoeff = [gaussCoeff, coeff];
        for peakNo = 1:length(coeff)
            
            f = coeff(strcat('peak',num2str(peakNo)));
            newX = slices(strcat('peak',num2str(peakNo)));
            newY = f(newX);
            plot(newX, newY, 'r+'); 
        end
        
        title(strcat("Curve fit for", foldernames(i).name));
        xlabel('pixels');
        ylabel('Amplitude');
        fileNo = fileNo + 1;
    catch ME
        switch ME.identifier
            case ''
                fprintf("not Valid directory\n")
            otherwise
                fprintf(ME.message)
                fprintf(ME.identifier)
                fprintf("\n")
                break
        end
    end

end

%Applying Gaussain Fit to ring peaks

function [values, sliceDictionary] = gaussFit(x, y, peaks)

fringeDictionary = containers.Map();
sliceDictionary = containers.Map();

% gaussFit = double.empty();

ringIndex = 1;
peakIndex = 2 ;
sliceIndex = 1;
while 1
    try

        ringEnds(ringIndex) = (peaks(peakIndex) + peaks(peakIndex+1)) / 2;

        
        xValues = x(sliceIndex:round(ringEnds(ringIndex)));
        yValues = y(sliceIndex:round(ringEnds(ringIndex)));
        
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
