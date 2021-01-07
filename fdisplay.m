function [] = fdisplay(folder, string, maxcol);

if isempty(folder); folder = pwd; end;

files = dir([folder '/bestfit*' string '*.txt']);

n = size(files);

nrows = ceil(n(1) / maxcol);

figure;

for j = 1:nrows
  for i = 1:maxcol

    c = (j - 1)*maxcol + i;

    m = load([folder '/' files(c).name]);

    subplot(nrows, maxcol, c), semilogy(m(:,1), m(:,2)); grid on;
    hold on; axis manual;
    semilogy(m(:,1), m(:,3), 'r', 'LineWidth', 2);

    hold off;

    if c == n(1); break; end;

  end
end
