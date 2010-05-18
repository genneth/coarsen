function [mu, p] = cdf_from_observations(file)

ms = load(file);
total = numel(ms);
average = sum(ms) / total;
counts = histc(ms, 0:max(ms));
p = cumsum(counts) / total;
mu = (0:max(ms)) / average;

end