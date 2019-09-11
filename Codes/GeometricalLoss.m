% Mon 13/09/2018
% Scr-F20
% geometrical loss
% Writer: Mojtaba Mansour Abadi
% Description: This file is used to calculae geometrical loss.


%% geometrical loss attenuation
a_src = Tx_Dir/norm(Tx_Dir);  % normalise the source direction vector
a_apr = Rx_Dir/norm(Rx_Dir);  % normalise the aperture direction vector
Pos_src_0 = Tx_Pos - Tx_Pos;  % move source position to the origin of the coordinate
Pos_apr_0 = Rx_Pos - Tx_Pos;  % move the aperture accordingly
% align, orient and translate the aperture based on source direction and
% orientation
V_link_0 = Pos_apr_0 - Pos_src_0;  % vector connecting source to receiver aqperture
Link_dist = norm(V_link_0);  % link distance (m)
[az, el, ~] = cart2sph(a_src(1), a_src(2), a_src(3));  % calculate the spherical coordinates; azimuth, elevation = 90 - theta, r
R_y = [[cos(-(pi/2 - el)),    0,     sin(-(pi/2 - el))];
       [0,                    1,     0                ];
       [-sin(-(pi/2 - el)),   0,     cos(-(pi/2 - el))]];  % rotate by '-(90 - el)' around y axis
R_z = [[cos(-az),    -sin(-az),     0];
       [sin(-az),     cos(-az),     0];
       [0,           0,             1]];  % rotate by '-az' around z axis;
R_dr = [[cosd(-Tx_Ori),    -sind(-Tx_Ori),  0];
       [sind(-Tx_Ori),    cosd(-Tx_Ori),   0];
       [0,              0,             1]];  % rotate by '-D_src' around z axis;
Apr_tr = R_dr*(R_y*(R_z*a_apr));  % rotatge the aperture direction vector to align z axis
a_apr_tr = Apr_tr/norm(Apr_tr);  % normalise the aligned aperture direction vector
Pos_apr_tr = R_dr*(R_y*(R_z*Pos_apr_0));  % transform the aperture position into new coordinate
V_link_tr = Pos_apr_tr - Pos_src_0;  % vector connecting source to receiver aqperture
a_link = V_link_tr/Link_dist;  % normalise the link vector
Aux_vec = dot(Apr_tr, a_link);  % dot product of two vectors
sai = 180 - acosd(Aux_vec);  % aperture angle relative to z axis (Deg)
% calculating the power
A_Rx_apr = pi*Rx_Ap_dia^2;  % area of the aperture
if(strcmp(TX_AP_Type, 'ANG'))  % angle of divergence is given
    theta_0_v = (theta_d_v*pi/180)/2;  % half vertical divergence angle (Rad)
    theta_0_h = (theta_d_h*pi/180)/2;  % half horizontal divergence angle (Rad)
    w_0_v = lambda/(pi*theta_0_v);  % aperture vertical radius (m)
    w_0_h = lambda/(pi*theta_0_h);  % aperture horizontal radius (m)
elseif(strcmp(TX_AP_Type, 'DIA'))  % aperture diameter is given
    w_0_v = w_tx_v/2;  % aperture vertical radius (m)
    w_0_h = w_tx_h/2;  % aperture horizontal radius (m)
else
    error('pick an aperture type for transmitter');  % print out an error message and exit
end
X = Pos_apr_tr(1);  % extract the x position of the receiver aperture (m)
Y = Pos_apr_tr(2);  % extract the y position of the receiver aperture (m)
Z = Pos_apr_tr(3);  % extract the z position of the receiver aperture (m)
z_0_v = pi*w_0_v^2/lambda;  % vertical Rayleigh range (m)
z_0_h = pi*w_0_h^2/lambda;  % horizontal Rayleigh range (m)
W_Rx_v = w_0_v*sqrt(1 + (Z/z_0_v)^2);  % vertical beam radius at receiver (m)
W_Rx_h = w_0_h*sqrt(1 + (Z/z_0_h)^2);  % horiontal beam radius at receiver (m)
I_uni = 1/(pi*W_Rx_v*W_Rx_h);  % uniform intensity
I_gus = 2/(pi*W_Rx_v*W_Rx_h)*exp(-2*X^2/W_Rx_v^2 -2*Y^2/W_Rx_h^2);  % Gaussian intensity
h_GL_uni = I_uni*A_Rx_apr*cosd(sai);  % geometrical loss; uniform intensity
h_GL_gus = I_gus*A_Rx_apr*cosd(sai);  % geometrical loss; Gaussian intensity
if(strcmp(Prop_Mod, 'UNI'))
    Geo_Loss_Sim = h_GL_uni;
elseif(strcmp(Prop_Mod, 'GUS'))
    Geo_Loss_Sim = h_GL_gus;
else
    error('pick an propagation model for beam');  % print out an error message and exit
end

if(abs(sai) > AFOV/2)  % check if the arrived beam direction is in the aperture FOV
    Geo_Loss_Sim = 0;
end

fprintf(['Geometrical loss:\tsai\t\t= %.2f Deg\n',...
    'Geometrical loss:\tLoss\t= %.2f dB\n'],...
    sai, 10*log10(Geo_Loss_Sim));  % print out the loss