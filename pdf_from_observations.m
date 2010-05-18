function [mu, p] = pdf_from_observations(file)

ms = load(file);
total = numel(ms);
average = sum(ms) / total;
counts = histc(ms, 1:max(ms));
p = counts / total * average;
mu = (1:max(ms)) / average;

end