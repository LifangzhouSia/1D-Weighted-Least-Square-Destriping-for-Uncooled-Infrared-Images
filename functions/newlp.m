function [final, down_bound, up_bound] = newlp(im)
% linear projection for 14bit image sequences
% input 14bit double image OR any bit
% output 8bit
im = double(im);
final = zeros(size(im));
p = size(im);
up_cut = 0.98;
down_cut = 0.02;
if length(p) == 2
    im_que = sort(im(:));
    up_bound = im_que(round(length(im_que)*up_cut));
    down_bound = im_que(round(length(im_que)*down_cut));
    im(im<=down_bound) = down_bound;
    im(im>=up_bound) = up_bound;
    dp = (im - min(im(:)))/(max(im(:))-min(im(:)));
    final = uint8(dp*255+0.5);
else
    for i = 1:p(3)
        imf = im(:,:,i);
        imf_que = sort(imf(:));
        up_bound = imf_que(round(length(imf_que)*up_cut));
        down_bound = imf_que(round(length(imf_que)*down_cut));
        imf(imf<=down_bound) = down_bound;
        imf(imf>=up_bound) = up_bound;
        dp = (imf - min(imf(:)))/(max(imf(:))-min(imf(:)));
        final(:,:,i) = uint8(dp*255+0.5);
    end
end
final = uint8(final);

end