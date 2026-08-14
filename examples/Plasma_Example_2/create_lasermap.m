clear; close all;

% constants
e = 1.602176565e-19;            % elementary charge (C)
m_e = 9.10938291e-31;           % electron mass (kg)
c = 2.99792458e8;		% speed of light (m/s)
E0 = 0.510998928e6;             % electron rest mass (eV/c)
eps_0 = 8.854187817e-12;        % vacuum permittivity (As/Vm)

filename = 'laser.dat';

w0 = 50e-6;                     % laser spot size in focus
a0 = 0.8;                       % laser normalized vector potential in focus
lambda_l = 800e-9;

zR = pi*w0^2/lambda_l;          % Rayleigh length

% build transverse vectors and construct mesh
% only +/- 1 rms of laser intensity is sampled,
% as we're only interested in the fields near the bunch
x_vec = linspace(-25e-6,25e-6,81);      
y_vec = linspace(-25e-6,25e-6,81);

plasma_start = 0;               % plasma start and end positions
plasma_end = 100e-3;

z_start_guiding = 10e-3;                       % laser focus position, here guiding starts
z_end_guiding = plasma_end - 10.e-3;           % coincides with end of plasma flat top region

z_vec = linspace(plasma_start,plasma_end,400);

[x,y,z] = meshgrid(x_vec,y_vec,z_vec);

% define laser envelope
% it is guided in plasma flat region, and follows free evolution in up and
% downramp
w = zeros(size(z)) + w0;       % set laser spot size to focal spot size everywhere
w(z<=z_start_guiding) = w0 * sqrt(1+(z(z<=z_start_guiding)-z_start_guiding).^2/zR^2);
w(z>=z_end_guiding) = w0 * sqrt(1+(z(z>=z_end_guiding)-z_end_guiding).^2/zR^2);

% calculate evolution of a^2 from evolution of spot size
asq = a0^2*w0^2./w.^2.*exp(-2*(x.^2+y.^2)./w.^2);

% reshape array to fit Astra 3D fieldmap requirements
asq_astra = reshape(permute(asq,[1 3 2]),length(y_vec)*length(z_vec),length(x_vec));

% construct header for Astra 3D fieldmap
header_x = [length(x_vec) x_vec(1) mean(diff(x_vec))];
header_y = [length(y_vec) y_vec(1) mean(diff(y_vec))];
header_z = [length(z_vec) z_vec(1) mean(diff(z_vec))];

save(filename,'header_x','header_y','header_z','asq_astra','-ascii','-double');
