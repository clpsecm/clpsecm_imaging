addpath(genpath('_algorithm/'));

d = DataProfile;
lines = read_clpsecm_data(d);

plot(lines.xticks,lines.current)