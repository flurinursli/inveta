function [] = space_search(fname, select);

% show the explored parameters space and misfit at every iteration. Use "select"
% to choose a specific parameter (leave empty to select all).
% e.g.: fdisplay('nasearch_2-4_MTI03.txt', []);


par = load(fname);

figure;

n = size(par);

if isempty(select)
  v = [1:n(2) - 1];
else
  v = select;
end

for i = 1:length(v)
  subplot(1, length(v), i), scatter(par(:, v(i)), par(:, n(2)), [], [1:n(1)], 'filled'); grid on;

  xlabel(['parameter #' num2str(i)])

  if i == 1; ylabel('Misfit'); end;

end

h = colorbar;
title(h, '# models explored');
