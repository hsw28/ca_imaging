function f = spikecenter(unsortedstructure, sortedstructure, dim)
  %finds brightest point and center point
  %dim should be the dimension you are expecting. will throw error if not correct
  %good are the cells manually selected


  extractedimages = unsortedstructure.cnmfeAnalysisOutput.extractedImages;

    if size(extractedimages,1)~=dim
      error('your dimension is not as expected');
    end

  good = sortedstructure.validCNMFE; %sorted cells
  temp = find(good==1);
  extractedimages = extractedimages(:, :, temp);


maxx = [];
maxy = [];
centerx = [];
centery = [];
for k=1:size(extractedimages,3)
    currentimage = extractedimages(:,:,k);
    currentimage = imgaussfilt(currentimage, 3); %smooth

    %find max values
    [x,y] = max(currentimage);
    maxx(end+1) = x;
    maxy(end+1) = y;

    %find center values
    [x,y] = find(currentimage>0);
    centerx(end+1) = mean(x);
    centery(end+1) = mean(y);

end

f = [maxx;maxy;centerx;centery];
