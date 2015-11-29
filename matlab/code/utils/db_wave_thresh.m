function im_rec = db_wave_thresh( im, lambda )

[im, r] = randshift(im);

dwtmode('ppd', 0);
wname = 'db4';

s = size(im);

im_rec = zeros(s);

for i = 1:prod(s(3:end));
    
    N = min( floor(log2( min(s(1:2) ) )), 3 );
    
    [coeff, book] = wavedec2( im(:,:,i), N, wname );
    coeff( prod(book(1,:))+1:end) = SoftThresh( coeff( prod(book(1,:))+1:end), lambda );
    
    
    im_rec(:,:,i) = waverec2( coeff, book, wname );
    
end

im_rec = randunshift( im_rec, r );