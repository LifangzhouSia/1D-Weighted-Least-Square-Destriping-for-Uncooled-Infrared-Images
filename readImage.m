function I = readImage(method)

pathstr = fileparts(which(method));
dirname = fullfile(pathstr, 'images', '*.png');
imglist = dir(dirname);

disp('  Available test images:');
disp(' ');
for k = 1:length(imglist)
  fprintf('  %d. %s\n', k, imglist(k).name);
end
disp(' ');

imnum = 0;
while (~isnumeric(imnum) || imnum<1 || imnum>length(imglist))
  imnum = input(sprintf('  Image to denoise (%d-%d): ', 1, length(imglist)), 's');
  imnum = sscanf(imnum, '%d');
end

imgname = fullfile(pathstr, 'images', imglist(imnum).name);


I = imread(imgname);

if length(size(I)) > 2
    I = double(rgb2gray(I));
else
    I = double(I);
end

end
