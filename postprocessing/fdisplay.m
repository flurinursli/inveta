function [] = fdisplay(folder, string, par);

% show observed (blue) and best fitting synthetic envelope (red). Plots can be
% arranged in "par(1)" rows and "par(2)" columns.
% e.g.: fdisplay('.', '8-16_MTI03_*', [3 4]);

nrows = par(1);
ncols = par(2);

if isempty(folder); folder = pwd; end;

files = dir([folder '/bestfit*' string '*.txt']);

n = size(files);

%nrows = ceil(n(1) / ncols);
nfigs = ceil(n(1) / (ncols * nrows));

for k = 1:nfigs

  figure;

  for j = 1:nrows
    for i = 1:ncols

      c = (k - 1)*nrows*ncols + (j - 1)*ncols + i;

      m = load([folder '/' files(c).name]);

      c = (j - 1)*ncols + i;

      subplot(nrows, ncols, c), semilogy(m(:,1), m(:,2)); grid on;
      hold on; axis manual;
      semilogy(m(:,1), m(:,3), 'r', 'LineWidth', 2);

      hold off;

      c = (k - 1)*nrows*ncols + (j - 1)*ncols + i;

      if c == n(1); break; end;

    end
  end

end
