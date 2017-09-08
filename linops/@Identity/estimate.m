function phase = estimate( op, im, r )

phase = angle( im * exp( 1i * r ) ) - r;
