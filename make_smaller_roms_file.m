function make_smaller_roms_file(input_filename, new_filename, xi, eta)
%MAKE_SMALLER_ROMS_FILE  subset history file for given xi and eta range
%  make_smaller_roms_file(input_filename, new_filename, xi, eta)
%
% variables: 'temp','salt','u_eastward','v_northward','zeta'
%
% - using constant z (not time-dependent vertical stretching) for now to
%   make it easier to integrate with cesium.
%
% Kevin Rosa
% June 7, 2020

fi = input_filename;
new_fi = new_filename;

%% determine dimension lengths
F = ncinfo(fi);
all_names = {F.Dimensions(:).Name};
all_lengths = [F.Dimensions(:).Length];

dimlen.time = all_lengths(strcmp(all_names, 'ocean_time'));
dimlen.N = all_lengths(strcmp(all_names, 'N'));

dimlen.xi = length(xi);
dimlen.eta = length(eta);

dimlen.one = 1;  % singleton dimension

%% create netcdf file and define dimensions
% Create/Open the file
ncid = netcdf.create(new_fi, 'CLOBBER');

% Define the dimensions
for dims = fieldnames(dimlen)'
    dim = dims{1};
    dimid.(dim) = netcdf.defDim(ncid, dim, dimlen.(dim));
end

%% define the variables
size_2d = [dimid.xi, dimid.eta];
size_3d = [dimid.xi, dimid.eta, dimid.time];
size_4d = [dimid.xi, dimid.eta, dimid.N, dimid.time];
size_z = [dimid.xi, dimid.eta, dimid.N];  % time-constant z

v_id.time = netcdf.defVar(ncid, 'ocean_time', 'double', [dimid.time]);

for vars = {'lon','lat','h'}
    var = vars{1};
    v_id.(var) = netcdf.defVar(ncid, var, 'double', size_2d);
end

for vars = {'zeta'}
    var = vars{1};
    v_id.(var) = netcdf.defVar(ncid, var, 'NC_FLOAT', size_3d);
end

for vars = {'u_eastward','v_northward','temp','salt'}
    var = vars{1};
    v_id.(var) = netcdf.defVar(ncid, var, 'NC_FLOAT', size_4d);
end

v_id.z = netcdf.defVar(ncid, 'z', 'double', size_z);


% Done defining the NetCdf
netcdf.endDef(ncid);


%% Read data and write to new netcdf
START = [xi(1) eta(1) 1 1];
COUNT = [length(xi) length(eta) Inf Inf];

ncwrite(new_fi, 'time', ncread(fi,'ocean_time'));

ncwrite(new_fi, 'lon', ncread(fi,'lon_rho',START(1:2),COUNT(1:2)));
ncwrite(new_fi, 'lat', ncread(fi,'lat_rho',START(1:2),COUNT(1:2)));
ncwrite(new_fi, 'h', ncread(fi,'h',START(1:2),COUNT(1:2)));

for vars = {'zeta'}
    var = vars{1};
    ncwrite(new_fi, var, ncread(fi,var,START([1,2,4]),COUNT([1,2,4])))
end

for vars = {'u_eastward','v_northward','temp','salt'}
    var = vars{1};
    ncwrite(new_fi, var, ncread(fi,var,START,COUNT))
end

%% compute z coordinates
for vars = {'Vtransform','Vstretching','theta_s','theta_b','Tcline','hc'}
    sig.(vars{1}) = ncread(fi, vars{1});
end
h = ncread(new_fi, 'h');
igrid = 1;  % rho points
z = set_depth(sig.Vtransform, sig.Vstretching, sig.theta_s, sig.theta_b, sig.hc, dimlen.N, ...
    igrid, h, 0, 0);

ncwrite(new_fi, 'z', z)

%% Metadata
% Re-enter define mode.
netcdf.reDef(ncid);

% Global attributes
global_varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid,global_varid,'creation_date',datestr(now));
netcdf.putAtt(ncid,global_varid,'created_by','Kevin Rosa');
netcdf.putAtt(ncid,global_varid,'source_file',fi);

% Variable attributes
for vars = fieldnames(v_id)'
    var = vars{1};

    i = find(strcmp({F.Variables(:).Name}, var));

    if strcmp(var, 'time')
        i = find(strcmp({F.Variables(:).Name}, 'ocean_time'));
    elseif strcmp(var, 'lon')
        i = find(strcmp({F.Variables(:).Name}, 'lon_rho'));
    elseif strcmp(var, 'lat_rho')
        i = find(strcmp({F.Variables(:).Name}, 'lat_rho'));
    elseif strcmp(var, 'mask_rho')
        i = find(strcmp({F.Variables(:).Name}, 'mask_rho'));
    end

    A = F.Variables(i).Attributes;
    for j = 1:length(A)
        netcdf.putAtt(ncid, v_id.(var), A(j).Name, A(j).Value)
    end

end
