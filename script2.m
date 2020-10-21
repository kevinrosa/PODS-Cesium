% Kevin Rosa
% Oct 21, 2020
addpath(genpath('~/repos/ROMS-Scripts'))

xi = 420:780;
eta = 550:1030;

t_inds = 1:365;
out_DIR = '/gpfs/data/epscor/dullman1/EPSCOR/2018_mod8_run1_NAM';
in_DIR = '/gpfs/data/epscor/krosa1/microplastics';

grid_fi = '/users/krosa1/projects/osom_2018/Input/osom_grid4_mindep_smlp_mod7.nc';

for t = t_inds
    fi_old = fullfile(out_DIR, sprintf('ocean_his_%04i.nc', t));
    time = nctime(fi_old,2,1);
    fi_new = fullfile(in_DIR, sprintf('osom_%s_xi%i-%i_eta%i-%i.nc',datestr(time,'yyyymmdd'),xi(1),xi(end),eta(1),eta(end)));
    
    make_smaller_roms_file_uv(fi_old, fi_new, grid_fi, xi, eta)
end
