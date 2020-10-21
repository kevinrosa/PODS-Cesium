function make_smaller_roms_file_uv(input_filename, new_filename, grid_fi, xi, eta)
%MAKE_SMALLER_ROMS_FILE  subset history file for given xi and eta range
%  make_smaller_roms_file(input_filename, new_filename, grid_fi, xi, eta)
%
% **
%
% Reworked Oct 13, 2020 to work with 3d velocities on u and v grids:
%   1. average u and v onto rho points.
%   2. use grid angle to rotate into u_eastward and v_northward.
%
% These runs also don't include grid info (lon, lat, angle, etc.) so the
% grid file will need to be provided separately.
%
% **
%
%
% variables: 'temp','salt','u_eastward','v_northward','zeta'
%
% - using constant z (not time-dependent vertical stretching) for now to
%   make it easier to integrate with cesium.
%
% --------------------------------------
% Converting to integers:
% 8bit:  +/-   128
% 16bit: +/- 32768
% 
%     Temp and Salt
% datatype: int16
% scale_factor: 1000
% add_offset: 15
% > range: -17 - +47
% > accuracy: 0.001
% 
%     Velocity                            Velocity (alternative)
% datatype: int16                     datatype: int8 
% scale_factor: 10000                 scale_factor: 67
% add_offset: 0                       add_offset: 0
% > range: -3.2 - +3.2                > range: -1.9 - +1.9
% > accuracy: 0.0001                  > accuracy: 0.015
% 
%     SSH (zeta)
% datatype: int8
% scale_factor: 50
% add_offset: 0
% > range: -2.5 - +2.5
% > accuracy: 0.02
% 
%     Longitude
% datatype: int16
% scale_factor: 20000
% add_offset: -71
% > range: -72.6 - -69.4
% > accuracy: 0.00005  (<5.6 m)
% 
%     Latitude
% datatype: int16
% scale_factor: 20000
% add_offset: 41
% > range: +39.4 - +42.6
% > accuracy: 0.00005  (5.6 m)
%
%     Depth (h) and Z
% datatype: int16
% scale_factor: 100
% add_offset: 0
% > range: -327 - +327
% > accuracy: 0.01  
%
%
% Info from Matlab ncwrite documentation:
%
% ncwrite function applies these attribute conventions in a sequence:
% 
%  1. Subtract the value of the add_offset attribute from vardata.
% 
%  2. Divide vardata by the value of the scale_factor attribute.
% 
%  3. Replace any NaN in vardata with the value contained in the _FillValue 
%     attribute. If this attribute does not exist, then ncwrite uses the 
%     fill value for this variable as specified by the NetCDF library.
% --------------------------------------
%
%
% Kevin Rosa
% June 7, 2020

%{ 


%}

fi = input_filename;
new_fi = new_filename;

F = ncinfo(fi);
G = ncinfo(grid_fi);

%% determine dimension lengths
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

v_id.time = netcdf.defVar(ncid, 'time', 'double', [dimid.time]);

for vars = {'lon','lat','h'}
    var = vars{1};
    v_id.(var) = netcdf.defVar(ncid, var, 'NC_SHORT', size_2d);
end

for vars = {'zeta'}
    var = vars{1};
    v_id.(var) = netcdf.defVar(ncid, var, 'NC_BYTE', size_3d);
end

for vars = {'u_eastward','v_northward','temp','salt'}
    var = vars{1};
    v_id.(var) = netcdf.defVar(ncid, var, 'NC_SHORT', size_4d);
end

v_id.z = netcdf.defVar(ncid, 'z', 'NC_SHORT', size_z);


%% Metadata

% Global attributes
global_varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid,global_varid,'creation_date',datestr(now));
netcdf.putAtt(ncid,global_varid,'created_by','Kevin Rosa');
netcdf.putAtt(ncid,global_varid,'source_file',fi);

% Variable attributes
for vars = fieldnames(v_id)'
    var = vars{1};
    
    switch var
        case 'time'
            i = find(strcmp({F.Variables(:).Name}, 'ocean_time'));
        case 'h'
            i = find(strcmp({G.Variables(:).Name}, 'h'));
        case 'lon'
            i = find(strcmp({G.Variables(:).Name}, 'lon_rho'));
        case 'lat'
            i = find(strcmp({G.Variables(:).Name}, 'lat_rho'));
        otherwise
            i = find(strcmp({F.Variables(:).Name}, var));
    end
    
    if ~isempty(i)
        switch var
            case {'h','lon','lat'}
                A = G.Variables(i).Attributes;
            otherwise
                A = F.Variables(i).Attributes;
        end
        
        for j = 1:length(A)
            netcdf.putAtt(ncid, v_id.(var), A(j).Name, A(j).Value)
        end
    end
end



%% scale factors, offsets, and fillvalues
for vars = {'temp','salt'}
    netcdf.putAtt(ncid, v_id.(vars{1}), 'scale_factor', 1/1000)
    netcdf.putAtt(ncid, v_id.(vars{1}), 'add_offset', 15)
end
for vars = {'u_eastward','v_northward'}
    netcdf.putAtt(ncid, v_id.(vars{1}), 'units', 'm/s')
    netcdf.putAtt(ncid, v_id.(vars{1}), 'scale_factor', 1/1000)
    netcdf.putAtt(ncid, v_id.(vars{1}), 'add_offset', 0)
    netcdf.putAtt(ncid, v_id.(vars{1}), 'units', 'm/s')
end
for vars = {'zeta'}
    netcdf.putAtt(ncid, v_id.(vars{1}), 'scale_factor', 1/50)
    netcdf.putAtt(ncid, v_id.(vars{1}), 'add_offset', 0)
end
for vars = {'lon'}
    netcdf.putAtt(ncid, v_id.(vars{1}), 'scale_factor', 1/20000)
    netcdf.putAtt(ncid, v_id.(vars{1}), 'add_offset', -71)
end
for vars = {'lat'}
    netcdf.putAtt(ncid, v_id.(vars{1}), 'scale_factor', 1/20000)
    netcdf.putAtt(ncid, v_id.(vars{1}), 'add_offset', 41)
end
for vars = {'h','z'}
    netcdf.putAtt(ncid, v_id.(vars{1}), 'units', 'meter')
    netcdf.putAtt(ncid, v_id.(vars{1}), 'scale_factor', 1/100)
    netcdf.putAtt(ncid, v_id.(vars{1}), 'add_offset', 0)
end


% Apply FillValues

% for int16 varibles:
for vars = {'temp','salt','u_eastward','v_northward','lon','lat','h','z'}
    netcdf.putAtt(ncid, v_id.(vars{1}), '_FillValue', int16(-32768))
end

% for int8 variables:
for vars = {'zeta'}
    netcdf.putAtt(ncid, v_id.(vars{1}), '_FillValue', int8(-128))
end

%% Exit define mode
netcdf.endDef(ncid);

%%
%
%
%% Read data and write to new netcdf
START = [xi(1) eta(1) 1 1];
COUNT = [length(xi) length(eta) Inf Inf];

ncwrite(new_fi, 'time', ncread(fi,'ocean_time'));

ncwrite(new_fi, 'lon', ncread(grid_fi,'lon_rho',START(1:2),COUNT(1:2)));
ncwrite(new_fi, 'lat', ncread(grid_fi,'lat_rho',START(1:2),COUNT(1:2)));
ncwrite(new_fi, 'h',   ncread(grid_fi,'h',START(1:2),COUNT(1:2)));

for vars = {'zeta'}
    var = vars{1};
    ncwrite(new_fi, var, ncread(fi,var,START([1,2,4]),COUNT([1,2,4])))
end

for vars = {'temp','salt'}
    var = vars{1};
    ncwrite(new_fi, var, ncread(fi,var,START,COUNT))
end

%% u and v: interpolate to rho-points
START_u = START;
COUNT_u = COUNT;    
START_u(1) = START(1)-1;
COUNT_u(1) = COUNT(1)+1;
u = ncread(fi,'u', START_u, COUNT_u);
u = (u(1:end-1,:,:,:) + u(2:end,:,:,:)) ./ 2;  % averaged onto rho points

START_v = START;
COUNT_v = COUNT;    
START_v(2) = START(2)-1;
COUNT_v(2) = COUNT(2)+1;
v = ncread(fi,'v', START_v, COUNT_v);
v = (v(:,1:end-1,:,:) + v(:,2:end,:,:)) ./ 2;

%% u and v: rotate 
ang = ncread(grid_fi,'angle',START(1:2),COUNT(1:2));  % angle between XI-axis and EAST (radians)

sz = size(u);
cos_ang = repmat(cos(ang), [1,1,sz(3),sz(4)]);  % 4D for multiplication
sin_ang = repmat(sin(ang), [1,1,sz(3),sz(4)]);  % 4D for multiplication

% dot products
u_eastward  = u .* cos_ang - v .* sin_ang;
v_northward = u .* sin_ang + v .* cos_ang;

%% u and v: write to netcdf
ncwrite(new_fi, 'u_eastward', u_eastward)
ncwrite(new_fi, 'v_northward', v_northward)


%% compute z coordinates
for vars = {'Vtransform','Vstretching','theta_s','theta_b','Tcline','hc'}
    sig.(vars{1}) = ncread(fi, vars{1});
end
h = ncread(new_fi, 'h');
igrid = 1;  % rho points
z = set_depth(sig.Vtransform, sig.Vstretching, sig.theta_s, sig.theta_b, sig.hc, dimlen.N, ...
    igrid, h, 0, 0);

ncwrite(new_fi, 'z', z)

end
