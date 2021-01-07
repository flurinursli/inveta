function [] = space_search(fname);

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
