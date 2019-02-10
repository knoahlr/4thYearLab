addpath("C:\Users\Noah Workstation\Desktop\P_PR_1\repo\attentionNN\lidar\atten\sim\noah_sim\");


%Opening and loading fringe data
fringeFolder = "D:\OneDrive - Carleton University\School\Courses Y4\Winter\PHYS 4007\Zeeman Effect\Data\Fringes_3";
DataFolder  = "..\Logs\Data\Fringes_3";
fitReport = "..\Logs\Data\Fringes_3\report.txt";
reportFile = fopen(fitReport, 'w');


foldernames = dir(DataFolder);
gaussCoeff = double.empty();
averageK = double.empty();
peakCenter = double.empty();
peakErrors = double.empty();
averageErrorK = double.empty();

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
        grid on;
        hold on;
        
        [coeff, slices, peakSpreads] = gaussFit(xValues, yValues, peaks);
%         averageV = [averageV, AnalysePeaks(peakCenter, foldernames(i).name)];
        
        fprintf(reportFile, "\n\n %s \n\n",foldernames(i).name);
        
        
        for peakNo = 1:length(coeff)

            equation = "f(x) = a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2)";
            
     
            f = coeff(strcat('peak',num2str(peakNo)));       
            coeffVals = coeffvalues(f);
            if coeffVals(2) > coeffVals(5)
                coeffVals = [coeffVals(4:6),coeffVals(1:3)];
            end
            
            Amplitude1 = num2str(coeffVals(1));
            center1 = num2str(coeffVals(2));
            error1 = num2str(coeffVals(3)/(2^0.5));
            temp = confint(f);
            centerConf1 = abs((temp(1,2) - temp(2,2)) /2);
            
            sigmaSpread1 = str2num(error1)/sqrt( peakSpreads(peakNo)/2);
            
            
            Amplitude2 = num2str(coeffVals(4));
            center2 = num2str(coeffVals(5));
            error2 = num2str(coeffVals(6)/(2^0.5));
            temp = confint(f);
            centerConf2 = abs((temp(1,4) - temp(2,4)) /2);
            
            sigmaSpread2 = str2num(error2)/sqrt(peakSpreads(peakNo)/2);
            
            if peakNo <= 4
                peakCenter = [peakCenter; str2num(center2),str2num(center1)];
                
                %peakErrors = [peakErrors; str2num(error2), str2num(error1)];
%                 peakErrors = [peakErrors; centerConf2, centerConf1];
                peakErrors = [peakErrors; sigmaSpread2, sigmaSpread1];
                
            else
                peakCenter = [peakCenter; str2num(center1),str2num(center2)];
                
                %peakErrors = [peakErrors; str2num(error1), str2num(error2)];
%                 peakErrors = [peakErrors; centerConf1, centerConf2];
                peakErrors = [peakErrors; sigmaSpread1, sigmaSpread2];
            end
            
                  
            fprintf(reportFile, "\n\n%s \n", strcat('Ring ',num2str(peakNo)));
            fprintf(reportFile, "Equation is %s\n", equation);
            fprintf(reportFile, "\nAmplitude1 = %s\nCenter1 = %s\nStandard Deviation = %s\nCenter1 Uncertainty = %s\n", Amplitude1, center1, error1, centerConf1);
            fprintf(reportFile, "\nAmplitude2 = %s\nCenter2 = %s\nStandard Deviation = %s\nCenter2 Uncertainty = %s\n", Amplitude2, center2, error2, centerConf2);
            
            newX = slices(strcat('peak',num2str(peakNo)));
            newY = f(newX);
            plot(newX, newY, 'r+'); 
        end
      
        
        title(strcat("Curve fit for", foldernames(i).name));
        xlabel('pixels');
        ylabel('Amplitude');
        fileNo = fileNo + 1;
        
        [deltaV, errorV] =  AnalysePeaks(peakCenter, peakErrors, foldernames(i).name);
        averageK = [averageK, deltaV];
        averageErrorK = [averageErrorK, errorV];
        peakCenter = double.empty();
        peakErrors = double.empty();
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
fclose(reportFile); 
%SLices data to isolate rings  and applys Gaussain Fit to ring peaks -
%gauss2






