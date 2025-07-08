function [emg_out, vel_out] = safeAlign(ratName, t_ref)

    [t_tmp, mean_emg, ~, mean_vel, ~] = alignEMGandVelocity(ratName);

    target_len = length(t_ref);
    cur_len = length(mean_emg);

    if cur_len == target_len
        emg_out = mean_emg;
        vel_out = mean_vel;
    elseif cur_len > target_len
        warning("%s trace too long; truncating", ratName);
        emg_out = mean_emg(1:target_len);
        vel_out = mean_vel(1:target_len);
    elseif cur_len < target_len
        warning("%s trace too short; padding", ratName);
        emg_out = [mean_emg, nan(1, target_len - cur_len)];
        vel_out = [mean_vel, nan(1, target_len - cur_len)];
    end
end
