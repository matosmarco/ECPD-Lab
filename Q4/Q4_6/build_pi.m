% Build Matrix pi
% function pi = build_pi(A,C,H)
% pi = zeros(1,H);
% for ii=1:H
%     pi(ii) = C*A^(ii);
% end
% pi = pi';
% end

function pi = build_pi(A, C, H)
% A dimensão do estado
n = size(A, 1);

% Boa prática: pré-alocar a matriz para o tamanho final H x n
pi = zeros(H, n);

for ii = 1:H
    % Guardamos o vetor linha na linha correspondente da matriz
    pi(ii, :) = C * (A^ii);
end
end