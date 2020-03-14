%%%%%%%%%%%% Authors: Javier Carnerero Cano    %%%%%%%%%%%%
%%%%%%%%%%%%          Vicente Gallardo Cabrera %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J = jacob(x_t, s)
    J = zeros(size(s, 2), size(s, 1));
    for i = 1 : size(J, 1)
        J(i, :) = ((s(:, i) - x_t) .* exp(-0.5 * norm(x_t - s(:, i))^2))';
    end
end