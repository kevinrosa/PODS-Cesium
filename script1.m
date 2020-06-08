% Kevin Rosa
% June 7, 2020
addpath(genpath('~/repos/ROMS-Scripts'))

xi = 420:780;
eta = 550:1030;

fi_old = '/users/krosa1/projects/osom_2018/Output/osom_2018_400_his_0010.nc';
fi_new = sprintf('%s_xi%i-%i_eta%i-%i_%s','/users/krosa1/projects/osom_2018/Output/osom_2018_400_his',xi(1),xi(end),eta(1),eta(end),'0010.nc');
make_smaller_roms_file(fi_old, fi_new, xi, eta)
