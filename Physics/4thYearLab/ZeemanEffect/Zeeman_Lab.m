
%for k = 0:13
for k = 10
    decimals = [2,7];
    whole = floor(48+k/2);
    index = strcat(string(whole), "p", string(decimals(mod(k,2)+1)));
    indRead = char(strcat(index, ".jpg"));
    sheet = char(index);
    pic = imread(indRead);
    %pic = imread('test2.jpg');
    newPic = pic([725:735],:,:);
    %newPic = pic(:,[850:900],:);
    %newPic = permute(newPic, [2,1,3]);
    imshow(newPic)
    a = sum(newPic);
    a = sum(a,3);
    aoriginal = a;
    a = movmean(a,10);
    aprime = diff(a);
    b = 1:length(a);
    d = 1:length(a)-1;
    c = [b',a'];
    
    check = 0;
    for n = 1:length(aprime)
        if aprime(n) < -50
            check = 1;
        end
        if check == 1 && aprime(n) > 0
            spot1 = n;
            break
        end
    end
    
    check = 0;
    for n = length(aprime):-1:1
        if aprime(n) > 50
            check = 1;
        end
        if check == 1 && aprime(n) < 0
            spot2 = n;
            break
        end
    end
    
    ashort = aoriginal(spot1:spot2);
    bshort = b(spot1:spot2);
    %figure
    %plot(a)
    %plot(d,aprime)
    
    figure
    %plot(ashort)
    plot(aoriginal)
    hold on
    %f = fit(bshort',ashort','gauss8');
    f = fit(b',aoriginal','gauss8');
    %vals = (f([spot1:spot2]));
    vals = (f([1:length(aoriginal)]));
    %plot(vals)
    coeffVals = coeffvalues(f);
    
    fringe_points = coeffVals([2:3:23]);
    fringe_errors = (coeffVals([3:3:24])/2).^0.5;
    
    check = 0;
    while check == 0
        check = 1;
        for n = 2:length(fringe_points)
            if fringe_points(n) < fringe_points(n-1)
                temp_point = fringe_points(n);
                fringe_points(n)= fringe_points(n-1);
                fringe_points(n-1) = temp_point;
                
                temp_point2 = fringe_errors(n);
                fringe_errors(n)= fringe_errors(n-1);
                fringe_errors(n-1) = temp_point2;
                
                check = 0;
            end
        end
    end
    fringe_points
    fringe_errors
    %xlswrite('C:\Users\Joseph\Documents\Carleton\Semester 8\PHYS 4007\Zeeman Lab\Zeeman Noah Joe\data.xlsx', c, sheet);
end