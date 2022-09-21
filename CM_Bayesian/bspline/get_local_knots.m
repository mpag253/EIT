function [local_knots, j] = get_local_knots(phi, global_knots)

    % get indices of global_knots points < phi
    % js = [p for p, q in enumerate(global_knots) if q < phi]
    lenm = length(global_knots);
    js = [];
    for p = 1:lenm, if global_knots(p) < phi, js = [js p]; end, end

    if isempty(js)
        j = 0;
        % local_knots = np.concatenate((global_knots[-2:] - 2*np.pi, global_knots[0:2]))
        local_knots = [global_knots(end-1:end)'-2*pi global_knots(1:2)'];

    else
        j = js(end);

        if j == 1
            % local_knots = np.concatenate(([global_knots[-1] - 2*np.pi], global_knots[0:3]))
            local_knots = [global_knots(end)-2*pi global_knots(1:3)'];

        elseif j == (lenm)
            % local_knots = np.concatenate((global_knots[-2:], global_knots[0:2] + 2*np.pi))
            local_knots = [global_knots(end-1:end)' global_knots(1:2)'+2*pi];

        elseif j == (lenm-1)
            % local_knots = np.concatenate((global_knots[-3:], [global_knots[0] + 2*np.pi]))
            local_knots = [global_knots(end-2:end)' global_knots(1)+2*pi];

        else
            % local_knots = global_knots[j - 1:j + 3]
            local_knots = global_knots(j-1:j+2);
        end
    end
end