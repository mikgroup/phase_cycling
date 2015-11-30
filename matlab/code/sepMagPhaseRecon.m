function [mag, phase] = sepMagPhaseRecon( ...
    ksp,...
    ESP_op, PFT_op, ...
    MAG_op, PHASE_op, ...
    MAG_thresh, PHASE_thresh, ...
    Mag_lambda, Phase_lambda, ...
    MAG_step, PHASE_step, ...
    nIter, do_phase_cycling, do_plot )

% Initialize print
if (~do_plot)
    it = 0;
    fprintf('Iteration %.3d/%.3d', it, nIter);
end

% Initialize mag and phase
resid = ESP_op' * (PFT_op' * ksp);
mag = abs( MAG_op' * resid );
Pp_init = angle( resid );
phase = PHASE_op' * Pp_init;


% Start iteration without regularization
for it = 1:nIter/5
    
    mag_old = mag;
    phase_old = phase;
    
    Mm = MAG_op * mag_old;
    expPp = exp(1j * (PHASE_op * phase_old));
    alpha_m = MAG_step;
    alpha_p = PHASE_step / (max(abs(Mm(:)))^2 + eps);
    
    % Update residue, mag, and phase
    resid = ESP_op' * (PFT_op' * ( ksp - PFT_op * (ESP_op * ( Mm .* expPp )) ));
    mag = mag_old + alpha_m * real( MAG_op' * ( conj(expPp) .* resid ) );
    phase = phase_old + alpha_p * imag( PHASE_op' * ( Mm .* conj(expPp) .* resid ) );
    
    % Plot
    if (do_plot)
        figure(31),
        clf
        imshow3( cat(2, abs(Mm) / max(abs(Mm(:))), abs(angle(expPp)) / max(abs(angle(expPp(:)))) ) );
        titlef(it);
        drawnow
    else
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bIteration %.3d/%.3d', it, nIter);
    end
end


% Start iteration
for it = (nIter/5+1):nIter
    
    
    mag_old = mag;
    phase_old = phase;
    
    Mm = MAG_op * mag_old;
    expPp = exp(1j * (PHASE_op * phase_old));
    alpha_m = MAG_step;
    alpha_p = PHASE_step / (max(abs(Mm(:)))^2 + eps);
    
    % Get random phase wraps
    if (do_phase_cycling)
        r = rand * 2*pi ;
        phase_wrap = PHASE_op' * (angle( exp(1j * (Pp_init + r) ) ) - Pp_init - r);
    else
        phase_wrap = zeros(size(phase));
    end
    
    % Update residue, mag, and phase
    resid = ESP_op' * (PFT_op' * ( ksp - PFT_op * (ESP_op * ( Mm .* expPp )) ));
    mag = MAG_thresh( mag_old + alpha_m * real( MAG_op' * ( conj(expPp) .* resid ) ), alpha_m * Mag_lambda );
    phase = PHASE_thresh( phase_old + alpha_p * imag( PHASE_op' * ( Mm .* conj(expPp) .* resid ) )...
        + phase_wrap, alpha_p * Phase_lambda ) - phase_wrap;
    
    
    % Plot
    if (do_plot)
        figure(31),
        clf
        imshow3( cat(2, abs(Mm) / max(abs(Mm(:))), abs(angle(expPp)) / max(abs(angle(expPp(:)))) ) );
        titlef(it);
        drawnow
    else
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bIteration %.3d/%.3d', it, nIter);
    end
end

if (~do_plot)
    fprintf('\n');
end