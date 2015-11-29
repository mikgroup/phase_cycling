function [mag, phase] = sepMagPhaseRecon( ...
    ksp, phase_init, ...
    ESP_op, PFT_op, ...
    MAG_op, PHASE_op, ...
    MAG_thresh, PHASE_thresh, ...
    Mag_lambda, Phase_lambda, ...
    MAG_step, PHASE_step, ...
    nIter, do_plot )

% Initialize print
if (~do_plot)
    it = 0;
    fprintf('Iteration %.3d/%.3d', it, nIter);
end

% Initialize mag and phase
resid = ESP_op' * (PFT_op' * ksp);
mag = abs( MAG_op \ resid ) * 0;
Pp_init = angle( resid .* conj(MAG_op * mag) );
phase = phase_init;

% Start iteration
for it = 1:nIter
    
    
    mag_old = mag;
    phase_old = phase;
    
    Mm = MAG_op * mag_old;
    expPp = exp(1j * (PHASE_op * phase_old));
    alpha_m = MAG_step;
    alpha_p = PHASE_step / (max(abs(Mm(:)))^2 + eps);
    
    % Get random phase wraps
%     r = rand * 2*pi ;
%     phase_wrap = angle( exp(1j * (Pp_init *2 + r) ) ) - 2*Pp_init - r;
    
    % Update residue, mag, and phase
    resid = ESP_op' * (PFT_op' * ( ksp - PFT_op * (ESP_op * ( Mm .* expPp )) ));
    mag = MAG_thresh( mag_old + alpha_m * real( MAG_op' * ( conj(expPp) .* resid ) ), Mag_lambda );
    phase = PHASE_thresh( phase_old + alpha_p * imag( PHASE_op' * ( Mm .* conj(expPp) .* resid ) ), Phase_lambda);
%     phase = PHASE_thresh( phase_old + alpha_p * imag( PHASE_op' * ( Mm .* conj(expPp) .* resid ) )...
%             + (PHASE_op' * phase_wrap), Phase_lambda ) - (PHASE_op' * phase_wrap);
    
        
    % Plot
    if (do_plot)
        figure(31),
        clf
        imshow3( cat(2,abs(mag)/max(abs(mag(:))),abs(phase)/pi)),
        titlef(it);
        drawnow
    else
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bIteration %.3d/%.3d', it, nIter);
    end
end

if (~do_plot)
    fprintf('\n');
end