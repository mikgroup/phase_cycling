function out = rand_phase(in,shift)

out = angle(exp(1i*(in+shift)));