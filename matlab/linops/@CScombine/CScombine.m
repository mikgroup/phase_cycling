function  res = CScombine( TE, FieldStrength, ppm)


res.adjoint = 0;
res.TE = TE;
res.FieldStrength = FieldStrength;
res.ppm = ppm;

ne = length(TE);
ns = length(ppm);

res.A = zeros(ne, ns);
GAMMA = 42.58; %in MHz/T

for s = 1:ns
    df = GAMMA * FieldStrength * ppm(s);
    for e = 1:ne
        res.A(e,s) = exp( -1i * 2*pi * TE(e) * df );
    end
end

res = class(res,'CScombine');

