function [OUT] = d1_WLS_Destriping(Insy, lambda, iter)

% initial image into [0, 1] range
initialFactor =  max(Insy(:));
L = Insy ./ initialFactor;
 
% initial the parameter
if(~exist('lambda', 'var')),
    lambda = 50;
end

if(~exist('iter', 'var')),
    iter = 1;
end

[hei, wid] = size(L);
u = L;
for t = 1:iter
    
    [Dire, miu] = edgeIndicator(u, 0.3);
    sigma_n = get_Sigma(u);

    for k = 1:hei
        u(k, :) = fGHS(u(k, :), ones(1, wid), lambda, 2*sigma_n);
    end
%     figure, imshow(newlp(u)), title('HF IMAGE')

    u = weightedLocalLinearFilter(u, L, Dire, floor(hei/6), 0.01);
        
end

% Output
OUT = u .* initialFactor;
s = (L - u) .* initialFactor;

end

%%%%%%%%%%%%%%% Horizontal Edge-preserving Smoothing Stage %%%%%%%%%%%%%%%

function u = fGHS(g, w, lambda, sigma_n)

[hei, wid] = size(g);
g = g(:);
w = w(:);

dg = diff(g);
dg = exp( - dg.^2 ./ (sigma_n^2) );
dga = - lambda*[0 dg.'].';
dgc = - lambda*[dg.' 0].';
dg  = w + lambda*([0 dg.'].' + [dg.' 0].');

u = linearInverseOperator(dga, dg, dgc, w.*g);

u = reshape(u, hei, wid);

end


function u = linearInverseOperator(a, b, c, f)
% using Gauss Elimination to inverse a 3-point laplacian matrix
DL = length(a);

% Calculate the forward computation
for k = 1:DL
    if k == 1
        c(k) = c(k)/b(k);
        f(k) = f(k)/b(k);
    else
        c(k) = c(k)/(b(k)- c(k-1)*a(k));
        f(k) = (f(k) - f(k-1)*a(k))/(b(k) - c(k-1)*a(k));
    end
end

% Calculate the backward computation
for k = DL:-1:1
    if k == DL
        u(k) = f(k);
    else
        u(k) = f(k) - c(k)*u(k+1);
    end
end
    
end

%%%%%%%%%%%%%%% Vertical Weighted Local Ridge Regression Stage %%%%%%%%%%%%%%%

function [Dire, m] = edgeIndicator(L, xi)

[hei, wid] = size(L);
r = 5;
Mean = (verticalBoxfilter(L.', r)).';
Mean2 = (verticalBoxfilter((L.*L).', r)).';
Var = Mean2 - Mean.*Mean;

m = mean(Var(:));
Dire = exp( - Var./(xi*m) ) + 1e-10;

end

function [u, a, b] = weightedLocalLinearFilter(x, y, w, r, eps)

[hei, wid] = size(x);

ww = verticalBoxfilter(w, r);
wx = verticalBoxfilter(w.*x, r)./ww;
wy = verticalBoxfilter(w.*y, r)./ww;
xwy = verticalBoxfilter(w.*x.*y, r)./ww;
wxx = verticalBoxfilter(w.*x.*x, r)./ww;
wyy = verticalBoxfilter(w.*y.*y, r)./ww;

a = (xwy - wx.*wy + eps) ./ (wxx - wx.*wx + eps);
b = wy - wx.*a;

mean_a =  verticalBoxfilter(a, r);
mean_b =  verticalBoxfilter(b, r);

u = (y - mean_b)./mean_a;

end

%%%%%%%%%%%%%%% Boxfilters and Some Functions %%%%%%%%%%%%%%%

function mI = verticalBoxfilter(I, r)

[hei, wid] = size(I);
Nv = D1boxfilter(ones(hei, 1), r);
mI = zeros(hei, wid);
for k = 1:wid
    mI(:, k) = D1boxfilter(I(:, k), r)./Nv;
end

end

function imDst = D1boxfilter(I, r)

imLine = I(:);
lth = length(imLine);
imDst = zeros(size(imLine));

imCum = cumsum(imLine, 1);
imDst(1:r+1) = imCum(1+r:2*r+1);
imDst(r+2:lth-r) = imCum(2*r+2:lth) - imCum(1:lth-2*r-1);
imDst(lth-r+1:lth) = repmat(imCum(lth), [r, 1]) - imCum(lth-2*r:lth-r-1);

end

function s = get_Sigma(I)

dx = diff(I, 1, 2);
s = 1.4826*median(abs(dx(:) - median(dx(:))));

end

