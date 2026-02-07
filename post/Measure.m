clear all;
format long;

cd(['saved_images']);
fid=fopen('R2.out','w');
for k=1:201
    fraction=zeros(6,1);
    image=imread(['poly_',num2str(k),'.jpg']);
    cropRect=[248,135,926,920];
    croppedImage=imcrop(image,cropRect);
    grayImage=rgb2gray(croppedImage);
    
    [counts,binValues]=imhist(grayImage);
    totalPixels=sum(counts);
    fractions=counts/totalPixels;
    [pks,locs]=findpeaks(fractions,'MinPeakHeight',0.0,'MinPeakDistance',5);
    if (k==1)
        [pks0,locs0]=findpeaks(fractions,'MinPeakHeight',0.01,'MinPeakDistance',5);
    end
    
    fraction=zeros(6,1);
    for i=1:6
        for j=1:size(pks,1)
            if (ismember(locs(j),locs0(i)-5:locs0(i)+5)&&pks(j)>0.0015)
                fraction(i)=sum(fractions(locs(j)-5:locs(j)+5));
            end
        end
    end
    fprintf(fid,'%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n',fraction(1),fraction(2),fraction(3),fraction(4),fraction(5),fraction(6));
end