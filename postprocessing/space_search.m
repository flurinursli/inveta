function [] = space_search(fname);

% show the explored parameters space and misfit at every iteration.
% e.g.: fdisplay('nasearch_2-4_MTI03.txt');


par = load(fname);

figure;

n = size(par);

for i = 1:n(2) - 1
  subplot(1, n(2) - 1, i), scatter(par(:, i), par(:, n(2)), [], [1:n(1)], 'filled'); grid on;

  xlabel(['parameter #' num2str(i)])

  if i == 1; ylabel('Misfit'); end;

end

h = colorbar;
title(h, '# models explored');
