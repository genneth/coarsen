function pure_voter_scaling

files = {'pure_voter-512-0.txt', ...
    'pure_voter-512-1.txt', ...
    'pure_voter-512-2.txt', ...
    'pure_voter-512-4.txt', ...
    'pure_voter-512-8.txt', ...
    'pure_voter-512-16.txt', ...
    'pure_voter-512-32.txt', ...
    'pure_voter-512-64.txt', ...
    'pure_voter-512-128.txt', ...
    'pure_voter-512-256.txt'
};

colours = [
 0.996078, 0.360784, 0.027451;
% 0.996078, 0.988235, 0.0352941;
 0.541176, 0.713725, 0.027451;
 0.145098, 0.435294, 0.384314;
 0.00784314, 0.509804, 0.929412;
 0.152941, 0.113725, 0.490196;
 0.470588, 0.262745, 0.584314;
 0.890196, 0.0117647, 0.490196;
 0.905882, 0.027451, 0.129412
];
ncolours = 8;

fh = figure;
gh = newplot(fh);
set(gh, 'NextPlot', 'add');
set(gh, 'XLim', [0 5], 'YScale', 'log', 'YLim', [10^-4 1]);
plot(gh, linspace(0, 5, 30), exp(-linspace(0,5,30)), '-k', 'LineWidth', 1.0);
av = zeros(size(files));
for i = 1:numel(files)
    [mu,p,av(i)] = cdf_from_observations(files{i});
    plot(gh, mu, 1-p, '-', 'Color', colours(mod(i-1,ncolours)+1,:) );
end

legend(gh, {'limit', '0','1','2','4','8','16', '32', '64', '128', '256'});

set(gh, 'FontName', 'Times', 'FontSize', 8);
xlabel(gh, 'rescaled clone size', 'FontName', 'Times', 'FontSize', 9);
ylabel(gh, 'cumulative distribution', 'FontName', 'Times', 'FontSize', 9);


%avh = newplot(figure);
%ts = [0 1 2 4 8 16 32 64];
%plot(avh, ts, pi() .* ts ./ log(ts), '-', ts, av, '+');
%set(avh, 'XLim', [0 64], 'YLim', [0 50]);

set(fh, 'PaperUnits', 'inches');
w = 5; h = 3.5;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

set(fh, 'Color', 'white');

drawnow;
print(fh, '-dpdf', 'pure_voter_scaling');

end
