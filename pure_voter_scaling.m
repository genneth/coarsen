function pure_voter_scaling

files = {'pure_voter-512-0.txt', ...
    'pure_voter-512-1.txt', ...
    'pure_voter-512-2.txt', ...
    'pure_voter-512-4.txt', ...
    'pure_voter-512-8.txt', ...
    'pure_voter-512-16.txt'};

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

gh = newplot(figure);
set(gh, 'NextPlot', 'add');
set(gh, 'XLim', [0 5], 'YScale', 'log', 'YLim', [10^-4 1]);
for i = 1:numel(files)
    [mu,p] = pdf_from_observations(files{i});
    plot(gh, mu, p, '-', 'Color', colours(i,:));
end

end
