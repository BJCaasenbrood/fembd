function [q, nu] = startState(s, x)
% Configuration and velocities of a state struct (its last step), or
% the initial configuration when x is empty.
    if isempty(x), q = s.q0;  nu = s.nu0;  return, end
    k = size(x.nu, 2);
    q = struct('r', x.r(:,:,k), 'R', x.R(:,:,:,k), 'e', x.e(:,k));  nu = x.nu(:,k);
end
