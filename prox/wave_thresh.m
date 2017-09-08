function P = wave_thresh(wname, N, lambda)
P = @(x, alpha) helper(x, alpha * lambda, wname, N);
end

function im_rec = helper(im, lambda, wname, N)

[im, r] = randshift(im);


dwtmode('ppd', 0);

s = size(im);

im_rec = zeros(s);

for i = 1:prod(s(3:end))
    
    
    [coeffr, bookr] = wavedec2( real(im(:,:,i)), N, wname );
    coeffr( prod(bookr(1,:))+1:end) = SoftThresh( coeffr( prod(bookr(1,:))+1:end), lambda );
    
    [coeffi, booki] = wavedec2( imag(im(:,:,i)), N, wname );
    coeffi( prod(booki(1,:))+1:end) = SoftThresh( coeffi( prod(booki(1,:))+1:end), lambda );

    
    im_rec(:,:,i) = waverec2( coeffr, bookr, wname ) + 1i * waverec2( coeffi, booki, wname );
    
end

im_rec = randunshift( im_rec, r );

end