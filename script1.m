% Kevin Rosa
% June 7, 2020
addpath(genpath('~/repos/ROMS-Scripts'))

xi = 420:780;
eta = 550:1030;

run_code = 'osom_2018_400';
t_inds = 10:40;
DIR = '/users/krosa1/projects/osom_2018/Output/';

for t = t_inds
    fi_old = fullfile(DIR, sprintf('%s_his_%04i.nc', run_code, t));
    time = nctime(fi_old,2,1);
    fi_new = fullfile(DIR, sprintf('osom_%s_xi%i-%i_eta%i-%i.nc',datestr(time,'yyyymmdd'),xi(1),xi(end),eta(1),eta(end)));
    make_smaller_roms_file(fi_old, fi_new, xi, eta)
end
