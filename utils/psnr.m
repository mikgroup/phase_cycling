function out = psnr(ref, compare)

out = 10 * log10(max(ref(:))^2 / mean((ref(:) - compare(:)).^2));