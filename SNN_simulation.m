function [ vars, records ] = SNN_simulation( vars, varargin )

% -------------------------------------------------------------------------

type_plasticity_off = find( strcmpi( varargin, 'plasticity_off' ) == 1 );

% type_triplet = find( strcmpi( varargin, 'triplet' ) == 1 );

type_records_on = find( strcmpi( varargin, 'records_on' ) == 1 );

% -------------------------------------------------------------------------

dt = vars.dt;
t_memory = vars.t_memory;
t_memory_L = t_memory / dt;

tau_E = vars.tau_E;
tau_I = vars.tau_I;
E_L_E = vars.E_L_E;
E_L_I = vars.E_L_I;
Delta_T_E = vars.Delta_T_E;
Delta_T_I = vars.Delta_T_I;
C = vars.C;
E_E = vars.E_E;
E_I = vars.E_I;
V_T = vars.V_T;
A_T = vars.A_T;
tau_T = vars.tau_T;
V_re = vars.V_re;
tau_abs = vars.tau_abs;
a_w = vars.a_w;
b_w = vars.b_w;
tau_w = vars.tau_w;
V_ap = vars.V_ap;
V_lb = vars.V_lb;

N_E = vars.N_E;
N_I = vars.N_I;
N_E_all = sum( N_E, 2 );
N_I_all = sum( N_I, 2 );
N_E_1 = N_E( 1 ) / 2;
N_E_2 = N_E( end ) / 2;
N_I_1 = N_I( 1 ) / 2;
N_I_2 = N_I( end ) / 2;
tau_r_E = vars.tau_r_E;
tau_d_E = vars.tau_d_E;
tau_r_I = vars.tau_r_I;
tau_d_I = vars.tau_d_I;
J_EE_min = vars.J_EE_min;
J_EE_max = vars.J_EE_max;
J_EI_min = vars.J_EI_min;
J_EI_max = vars.J_EI_max;
J_ext_EE = vars.J_ext_EE;
J_ext_IE = vars.J_ext_IE;
period_scaling = vars.period_scaling;

A_LTD = vars.A_LTD;
A_LTP = vars.A_LTP;
theta_LTD = vars.theta_LTD;
theta_LTP = vars.theta_LTP;
tau_u = vars.tau_u;
tau_v = vars.tau_v;
tau_x = vars.tau_x;

tau_y = vars.tau_y;
eta = vars.eta;
r_O = vars.r_O / 1000;

% triplet_tau_plus = vars.triplet_tau_plus;
% triplet_tau_minus = vars.triplet_tau_minus;
% triplet_tau_x = vars.triplet_tau_x;
% triplet_tau_y = vars.triplet_tau_y;
% triplet_A_2_plus = vars.triplet_A_2_plus;
% triplet_A_3_plus = vars.triplet_A_3_plus;
% triplet_A_2_minus = vars.triplet_A_2_minus;
% triplet_A_3_minus = vars.triplet_A_3_minus;

% -------------------------------------------------------------------------

idx_J_EE = vars.idx_J_EE;
idx_J_EI = vars.idx_J_EI;
idx_J_IE = vars.idx_J_IE;
idx_J_II = vars.idx_J_II;
% idx_J_GII = vars.idx_J_GII;

idx_J_EE_1 = vars.idx_J_EE_1;
idx_J_EI_1 = vars.idx_J_EI_1;
idx_J_IE_1 = vars.idx_J_IE_1;
idx_J_II_1 = vars.idx_J_II_1;

J_EE = vars.J_EE;
J_EI = vars.J_EI;
J_IE = vars.J_IE;
J_II = vars.J_II;
% J_GII = vars.J_GII;

init_J_EE = vars.init_J_EE;
connects = vars.connects;
idxMat_E = vars.idxMat_E;

% -------------------------------------------------------------------------

if vars.mode == 1
    J_EE_1 = J_EE;
    J_EI_1 = J_EI;
    J_IE_1 = J_IE;
    J_II_1 = J_II;

    J_EE_1( ~idx_J_EE_1 ) = 0;
    J_EI_1( ~idx_J_EI_1 ) = 0;
    J_IE_1( ~idx_J_IE_1 ) = 0;
    J_II_1( ~idx_J_II_1 ) = 0;

    ratio_EE = sum( idx_J_EE_1( : ), 1 ) / sum( idx_J_EE( : ), 1 );
    ratio_EI = sum( idx_J_EI_1( : ), 1 ) / sum( idx_J_EI( : ), 1 );
    ratio_IE = sum( idx_J_IE_1( : ), 1 ) / sum( idx_J_IE( : ), 1 );
    ratio_II = sum( idx_J_II_1( : ), 1 ) / sum( idx_J_II( : ), 1 );
    
    J_EE_1 = J_EE_1 / ratio_EE;
    J_EI_1 = J_EI_1 / ratio_EI;
    J_IE_1 = J_IE_1 / ratio_IE;
    J_II_1 = J_II_1 / ratio_II;

    J_EE_min = J_EE_min / ratio_EE;
    J_EE_max = J_EE_max / ratio_EE;
    J_EI_min = J_EI_min / ratio_EI;
    J_EI_max = J_EI_max / ratio_EI;

    J_EE_O = vars.J_EE_O;
    J_EE_O_1 = J_EE_O / ratio_EE;
    init_J_EE_1 = J_EE_O_1 * idx_J_EE_1;
end

% -------------------------------------------------------------------------

inputs = vars.inputs;
ext_stim_E = vars.ext_stim_E;
ext_stim_I = vars.ext_stim_I;

t_domain = 1 : size( ext_stim_E, 2 );
t_domain_scaling = round( period_scaling / dt : period_scaling / dt : length( t_domain ) );
t_domain_scaling = ismember( t_domain, t_domain_scaling );
t_domain_display = round( linspace( 0, length( t_domain ), 6 ) );
t_domain_display = t_domain_display( 2 : end );
t_domain_display = ismember( t_domain, t_domain_display );

% -------------------------------------------------------------------------

V_E = vars.V_E;
V_I = vars.V_I;
Vu_E = vars.Vu_E;
Vv_E = vars.Vv_E;

s_E = vars.s_E;
sx_E = vars.sx_E;
sy_E = vars.sy_E;
s_I = vars.s_I;
sy_I = vars.sy_I;

V_T_E = vars.V_T_E;
ww_E = vars.ww_E;

s_ext_EE = vars.s_ext_EE;
s_ext_IE = vars.s_ext_IE;

ct_abs_E = vars.ct_abs_E;
ct_abs_I = vars.ct_abs_I;

% triplet_r_1 = vars.triplet_r_1;
% triplet_r_2 = vars.triplet_r_2;
% triplet_o_1 = vars.triplet_o_1;
% triplet_o_2 = vars.triplet_o_2;

% -------------------------------------------------------------------------

if isempty( V_E )
    V_E = E_L_E * ones( N_E_all, t_memory_L );
end
if isempty( V_I )
    V_I = E_L_I * ones( N_I_all, t_memory_L );
end
if isempty( Vu_E )
    Vu_E = E_L_E * ones( N_E_all, t_memory_L );
end
if isempty( Vv_E )
    Vv_E = E_L_E * ones( N_E_all, t_memory_L );
end

if isempty( s_E )
    s_E = zeros( N_E_all, t_memory_L );
end
if isempty( sx_E )
    sx_E = zeros( N_E_all, t_memory_L );
end
if isempty( sy_E )
    sy_E = zeros( N_E_all, t_memory_L );
end
if isempty( s_I )
    s_I = zeros( N_I_all, t_memory_L );
end
if isempty( sy_I )
    sy_I = zeros( N_I_all, t_memory_L );
end

if isempty( V_T_E )
    V_T_E = V_T * ones( N_E_all, 1 );
end
if isempty( ww_E )
    ww_E = 0 * ones( N_E_all, 1 );
end

if isempty( s_ext_EE )
    s_ext_EE = zeros( N_E_all, t_memory_L );
end
if isempty( s_ext_IE )
    s_ext_IE = zeros( N_I_all, t_memory_L );
end

if isempty( ct_abs_E )
    ct_abs_E = zeros( N_E_all, 1 );
end
if isempty( ct_abs_I )
    ct_abs_I = zeros( N_I_all, 1 );
end

% if isempty( triplet_r_1 )
%     triplet_r_1 = zeros( N_E_all, t_memory_L );
% end
% if isempty( triplet_r_2 )
%     triplet_r_2 = zeros( N_E_all, t_memory_L );
% end
% if isempty( triplet_o_1 )
%     triplet_o_1 = zeros( N_E_all, t_memory_L );
% end
% if isempty( triplet_o_2 )
%     triplet_o_2 = zeros( N_E_all, t_memory_L );
% end

% -------------------------------------------------------------------------

t_domain_memory = 0 : dt : t_memory - dt;
F_E = ( 1 / ( tau_d_E - tau_r_E ) ) * ( exp( -t_domain_memory / tau_d_E ) - exp( -t_domain_memory / tau_r_E ) );
F_I = ( 1 / ( tau_d_I - tau_r_I ) ) * ( exp( -t_domain_memory / tau_d_I ) - exp( -t_domain_memory / tau_r_I ) );
t_domain_memory = fliplr( -t_domain_memory );
F_E = fliplr( F_E );
F_I = fliplr( F_I );
% figure
% hold on
% plot( t_domain_memory, F_E, 'r' )
% plot( t_domain_memory, F_I, 'b' )

% -------------------------------------------------------------------------

temp_V_E = V_E( :, end );
temp_V_I = V_I( :, end );
temp_Vu_E = Vu_E( :, end );
temp_Vv_E = Vv_E( :, end );

temp_s_E = s_E( :, end );
temp_sx_E = sx_E( :, end );
temp_sy_E = sy_E( :, end );
temp_s_I = s_I( :, end );
temp_sy_I = sy_I( :, end );

temp_s_ext_EE = s_ext_EE( :, end );
temp_s_ext_IE = s_ext_IE( :, end );

temp_R_Vu_E = zeros( N_E_all, 1 );
temp_R_V_E = zeros( N_E_all, 1 );
temp_R_Vv_E = zeros( N_E_all, 1 );

g_EE = zeros( N_E_all, 1 );
g_EI = zeros( N_E_all, 1 );
g_IE = zeros( N_I_all, 1 );
g_II = zeros( N_I_all, 1 );
% g_GII = zeros( N_I_all, 1 );

% temp_triplet_r_1 = triplet_r_1( :, end );
% temp_triplet_r_2 = triplet_r_2( :, end );
% temp_triplet_o_1 = triplet_o_1( :, end );
% temp_triplet_o_2 = triplet_o_2( :, end );

% temp_triplet_pre = zeros( N_E_all, N_E_all );
% temp_triplet_post = zeros( N_E_all, N_E_all );

% -------------------------------------------------------------------------

if ~isempty( type_records_on )
    clear records
    records.s_E = false( N_E_all, length( t_domain ) );
    records.s_I = false( N_I_all, length( t_domain ) );
end

% -------------------------------------------------------------------------

ct_domain_display = 0;
for t = t_domain

    % -------------------------------------------------------------------------
    % voltages

    temp_s_ext_EE( :, 1 ) = ext_stim_E( :, t );
    temp_s_ext_IE( :, 1 ) = ext_stim_I( :, t );
    s_ext_EE( :, 1 : end - 1 ) = s_ext_EE( :, 2 : end );
    s_ext_IE( :, 1 : end - 1 ) = s_ext_IE( :, 2 : end );
    s_ext_EE( :, end ) = temp_s_ext_EE;
    s_ext_IE( :, end ) = temp_s_ext_IE;

    if vars.mode < 1
        g_EE = sum( F_E .* ( ( J_ext_EE * s_ext_EE ) + ( J_EE * s_E ) ), 2 );
        g_EI = sum( F_I .* ( J_EI * s_I ), 2 );
        g_IE = sum( F_E .* ( ( J_ext_IE * s_ext_IE ) + ( J_IE * s_E ) ), 2 );
        g_II = sum( F_I .* ( J_II * s_I ), 2 );
        % g_GII = sum( J_GII .* ( transpose( temp_V_I ) - temp_V_I ), 2 );
    elseif vars.mode == 1
        g_EE = sum( F_E .* ( ( J_ext_EE * s_ext_EE ) + ( J_EE_1 * s_E ) ), 2 );
        g_EI = sum( F_I .* ( J_EI_1 * s_I ), 2 );
        g_IE = sum( F_E .* ( ( J_ext_IE * s_ext_IE ) + ( J_IE_1 * s_E ) ), 2 );
        g_II = sum( F_I .* ( J_II_1 * s_I ), 2 );
    end

    V_T_E = V_T_E + dt * ( ( 1 / tau_T ) * ( V_T - V_T_E ) );
    ww_E = ww_E + dt * ( ( 1 / tau_w ) * ( ( a_w * ( temp_V_E - E_L_E ) ) - ww_E ) );
    ww_E = ww_E + ( b_w * temp_s_E );

    temp_V_E = temp_V_E + dt * ( ( 1 / tau_E ) * ( E_L_E - temp_V_E + Delta_T_E * exp( ( temp_V_E - V_T_E ) / Delta_T_E ) ) + ( g_EE / C ) .* ( E_E - temp_V_E ) + ( g_EI / C ) .* ( E_I - temp_V_E ) - ( ww_E / C ) );
    % temp_V_I = temp_V_I + dt * ( ( 1 / tau_I ) * ( E_L_I - temp_V_I + Delta_T_I * exp( ( temp_V_I - V_T ) / Delta_T_I ) ) + ( g_IE / C ) .* ( E_E - temp_V_I ) + ( g_II / C ) .* ( E_I - temp_V_I ) + ( g_GII / C ) );
    temp_V_I = temp_V_I + dt * ( ( 1 / tau_I ) * ( E_L_I - temp_V_I + Delta_T_I * exp( ( temp_V_I - V_T ) / Delta_T_I ) ) + ( g_IE / C ) .* ( E_E - temp_V_I ) + ( g_II / C ) .* ( E_I - temp_V_I ) );

    % -------------------------------------------------------------------------
    % spikes & re-voltages

    ct_abs_E( ct_abs_E > 0 ) = ct_abs_E( ct_abs_E > 0 ) - dt;
    ct_abs_I( ct_abs_I > 0 ) = ct_abs_I( ct_abs_I > 0 ) - dt;

    temp_V_E( temp_V_E <= V_lb ) = V_lb;
    temp_V_I( temp_V_I <= V_lb ) = V_lb;

    temp_V_E( temp_V_E >= V_ap ) = V_ap;
    temp_V_I( temp_V_I >= V_ap ) = V_ap;

    temp_s_E = ct_abs_E <= 0 & temp_V_E >= V_ap;
    temp_s_I = ct_abs_I <= 0 & temp_V_I >= V_ap;

    ct_abs_E( temp_s_E ) = tau_abs;
    ct_abs_I( temp_s_I ) = tau_abs;

    temp_V_E( ct_abs_E > 0 ) = V_re;
    temp_V_I( ct_abs_I > 0 ) = V_re;

    V_E( :, 1 : end - 1 ) = V_E( :, 2 : end );
    V_I( :, 1 : end - 1 ) = V_I( :, 2 : end );
    s_E( :, 1 : end - 1 ) = s_E( :, 2 : end );
    s_I( :, 1 : end - 1 ) = s_I( :, 2 : end );
    V_E( :, end ) = temp_V_E;
    V_I( :, end ) = temp_V_I;
    s_E( :, end ) = temp_s_E;
    s_I( :, end ) = temp_s_I;

    V_T_E( temp_s_E ) = V_T + A_T;

    if ~isempty( type_records_on )
        records.s_E( :, t ) = temp_s_E;
        records.s_I( :, t ) = temp_s_I;
    end

    % -------------------------------------------------------------------------
    % low-pass filtering

    temp_Vu_E = temp_Vu_E + ( dt / ( dt + tau_u ) ) * ( temp_V_E - temp_Vu_E );
    temp_Vv_E = temp_Vv_E + ( dt / ( dt + tau_v ) ) * ( temp_V_E - temp_Vv_E );
    temp_sx_E = temp_sx_E + ( dt / ( dt + tau_x ) ) * ( temp_s_E - temp_sx_E );
    temp_sy_E = temp_sy_E + ( dt / ( dt + tau_y ) ) * ( temp_s_E - temp_sy_E );
    temp_sy_I = temp_sy_I + ( dt / ( dt + tau_y ) ) * ( temp_s_I - temp_sy_I );
    Vu_E( :, 1 : end - 1 ) = Vu_E( :, 2 : end );
    Vv_E( :, 1 : end - 1 ) = Vv_E( :, 2 : end );
    sx_E( :, 1 : end - 1 ) = sx_E( :, 2 : end );
    sy_E( :, 1 : end - 1 ) = sy_E( :, 2 : end );
    sy_I( :, 1 : end - 1 ) = sy_I( :, 2 : end );
    Vu_E( :, end ) = temp_Vu_E;
    Vv_E( :, end ) = temp_Vv_E;
    sx_E( :, end ) = temp_sx_E;
    sy_E( :, end ) = temp_sy_E;
    sy_I( :, end ) = temp_sy_I;

    % -------------------------------------------------------------------------
    % triplet variable update

    %     temp_triplet_r_1 = temp_triplet_r_1 + dt * ( ( 1 / triplet_tau_plus ) * ( -1 ) * temp_triplet_r_1 );
    %     temp_triplet_r_2 = temp_triplet_r_2 + dt * ( ( 1 / triplet_tau_x ) * ( -1 ) * temp_triplet_r_2 );
    %     temp_triplet_o_1 = temp_triplet_o_1 + dt * ( ( 1 / triplet_tau_minus ) * ( -1 ) * temp_triplet_o_1 );
    %     temp_triplet_o_2 = temp_triplet_o_2 + dt * ( ( 1 / triplet_tau_y ) * ( -1 ) * temp_triplet_o_2 );
    %     temp_triplet_r_1( temp_s_E ) = temp_triplet_r_1( temp_s_E ) + 1;
    %     temp_triplet_r_2( temp_s_E ) = temp_triplet_r_2( temp_s_E ) + 1;
    %     temp_triplet_o_1( temp_s_E ) = temp_triplet_o_1( temp_s_E ) + 1;
    %     temp_triplet_o_2( temp_s_E ) = temp_triplet_o_2( temp_s_E ) + 1;
    %     triplet_r_1( :, 1 : end - 1 ) = triplet_r_1( :, 2 : end );
    %     triplet_r_2( :, 1 : end - 1 ) = triplet_r_2( :, 2 : end );
    %     triplet_o_1( :, 1 : end - 1 ) = triplet_o_1( :, 2 : end );
    %     triplet_o_2( :, 1 : end - 1 ) = triplet_o_2( :, 2 : end );
    %     triplet_r_1( :, end ) = temp_triplet_r_1;
    %     triplet_r_2( :, end ) = temp_triplet_r_2;
    %     triplet_o_1( :, end ) = temp_triplet_o_1;
    %     triplet_o_2( :, end ) = temp_triplet_o_2;

    % -------------------------------------------------------------------------
    % synaptic plasticity

    if isempty( type_plasticity_off )
        if vars.mode < 1
            % if isempty( type_triplet )% voltage-based plasticity
            temp_R_Vu_E = temp_Vu_E - theta_LTD;
            temp_R_V_E = temp_V_E - theta_LTP;
            temp_R_Vv_E = temp_Vv_E - theta_LTD;
            temp_R_Vu_E( temp_R_Vu_E < 0 ) = 0;
            temp_R_V_E( temp_R_V_E < 0 ) = 0;
            temp_R_Vv_E( temp_R_Vv_E < 0 ) = 0;
            J_EE = J_EE + dt * idx_J_EE .* ( - A_LTD * temp_R_Vu_E * transpose( temp_s_E ) + A_LTP * ( temp_R_V_E .* temp_R_Vv_E ) * transpose( temp_sx_E ) );
            J_EE( idx_J_EE & J_EE < J_EE_min ) = J_EE_min;
            J_EE( idx_J_EE & J_EE > J_EE_max ) = J_EE_max;
            J_EE_1( idx_J_EE_1 ) = J_EE( idx_J_EE_1 );
            % elseif ~isempty( type_triplet )% triplet plasticity
            %             temp_triplet_pre = temp_triplet_o_1 .* ( triplet_A_2_minus + triplet_A_3_minus .* transpose( temp_triplet_r_2 ) );
            %             temp_triplet_post = transpose( temp_triplet_r_1 ) .* ( triplet_A_2_plus + triplet_A_3_plus .* temp_triplet_o_2 );
            %             J_EE( :, transpose( temp_s_E ) ) = J_EE( :, transpose( temp_s_E ) ) - temp_triplet_pre( :, transpose( temp_s_E ) );
            %             J_EE( temp_s_E, : ) = J_EE( temp_s_E, : ) + temp_triplet_post( temp_s_E, : );
            % end
            %             J_EI = J_EI + eta * idx_J_EI .* ( ( temp_sy_E - 2 * r_O * tau_y ) * transpose( temp_s_I ) );
            %             J_EI = J_EI + eta * idx_J_EI .* ( temp_s_E * transpose( temp_sy_I ) );
            %             J_EI( idx_J_EI & J_EI < J_EI_min ) = J_EI_min;
            %             J_EI( idx_J_EI & J_EI > J_EI_max ) = J_EI_max;
        elseif vars.mode == 1
            % if isempty( type_triplet )% voltage-based plasticity
            temp_R_Vu_E = temp_Vu_E - theta_LTD;
            temp_R_V_E = temp_V_E - theta_LTP;
            temp_R_Vv_E = temp_Vv_E - theta_LTD;
            temp_R_Vu_E( temp_R_Vu_E < 0 ) = 0;
            temp_R_V_E( temp_R_V_E < 0 ) = 0;
            temp_R_Vv_E( temp_R_Vv_E < 0 ) = 0;
            J_EE_1 = J_EE_1 + dt * idx_J_EE_1 .* ( - A_LTD * temp_R_Vu_E * transpose( temp_s_E ) + A_LTP * ( temp_R_V_E .* temp_R_Vv_E ) * transpose( temp_sx_E ) );
            J_EE_1( idx_J_EE_1 & J_EE_1 < J_EE_min ) = J_EE_min;
            J_EE_1( idx_J_EE_1 & J_EE_1 > J_EE_max ) = J_EE_max;
            J_EE( idx_J_EE_1 ) = J_EE_1( idx_J_EE_1 );
            % elseif ~isempty( type_triplet )% triplet plasticity
            %             temp_triplet_pre = temp_triplet_o_1 .* ( triplet_A_2_minus + triplet_A_3_minus .* transpose( temp_triplet_r_2 ) );
            %             temp_triplet_post = transpose( temp_triplet_r_1 ) .* ( triplet_A_2_plus + triplet_A_3_plus .* temp_triplet_o_2 );
            %             J_EE( :, transpose( temp_s_E ) ) = J_EE( :, transpose( temp_s_E ) ) - temp_triplet_pre( :, transpose( temp_s_E ) );
            %             J_EE( temp_s_E, : ) = J_EE( temp_s_E, : ) + temp_triplet_post( temp_s_E, : );
            % end
            %             J_EI_1 = J_EI_1 + eta * idx_J_EI_1 .* ( ( temp_sy_E - 2 * r_O * tau_y ) * transpose( temp_s_I ) );
            %             J_EI_1 = J_EI_1 + eta * idx_J_EI_1 .* ( temp_s_E * transpose( temp_sy_I ) );
            %             J_EI_1( idx_J_EI_1 & J_EI_1 < J_EI_min ) = J_EI_min;
            %             J_EI_1( idx_J_EI_1 & J_EI_1 > J_EI_max ) = J_EI_max;
        end

        if t_domain_scaling( t )
            if vars.mode < 1
                for h1 = 1 : length( N_E )
                    for h2 = 1 : length( N_E )
                        if vars.connects( h1, h2 ) > 0
                            % J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) = J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) - idx_J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) .* ( sum( J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) - init_J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ), 2 ) ./ sum( idx_J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ), 2 ) );
                            J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) = vars.connects( h1, h2 ) * ( J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) - idx_J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) .* ( sum( J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) - init_J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ), 2 ) ./ sum( idx_J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ), 2 ) ) );
                        end
                    end
                end
                J_EE( isnan( J_EE ) ) = 0;
                J_EE( idx_J_EE & J_EE < J_EE_min ) = J_EE_min;
                J_EE( idx_J_EE & J_EE > J_EE_max ) = J_EE_max;
            elseif vars.mode == 1
                for h1 = 1 : length( N_E )
                    for h2 = 1 : length( N_E )
                        if vars.connects( h1, h2 ) > 0 && h1 ~= h2
                            J_EE_1( idxMat_E( h1, : ), idxMat_E( h2, : ) ) = J_EE_1( idxMat_E( h1, : ), idxMat_E( h2, : ) ) - idx_J_EE_1( idxMat_E( h1, : ), idxMat_E( h2, : ) ) .* ( sum( J_EE_1( idxMat_E( h1, : ), idxMat_E( h2, : ) ) - init_J_EE_1( idxMat_E( h1, : ), idxMat_E( h2, : ) ), 2 ) ./ sum( idx_J_EE_1( idxMat_E( h1, : ), idxMat_E( h2, : ) ), 2 ) );
                        end
                    end
                end
                J_EE_1( isnan( J_EE_1 ) ) = 0;
                J_EE_1( idx_J_EE_1 & J_EE_1 < J_EE_min ) = J_EE_min;
                J_EE_1( idx_J_EE_1 & J_EE_1 > J_EE_max ) = J_EE_max;
            end
        end
    end

    % -------------------------------------------------------------------------
    % progress display

    %     if t_domain_display( t )
    %         ct_domain_display = ct_domain_display + 20;
    %         disp( [ num2str( ct_domain_display ), '% complete!' ] )
    %     end

    % -------------------------------------------------------------------------

end

% -------------------------------------------------------------------------

if vars.mode == 1
    J_EE_1 = J_EE_1 * ratio_EE;
    J_EI_1 = J_EI_1 * ratio_EI;

    J_EE_min = J_EE_min * ratio_EE;
    J_EE_max = J_EE_max * ratio_EE;
    J_EI_min = J_EI_min * ratio_EI;
    J_EI_max = J_EI_max * ratio_EI;

    J_EE( idx_J_EE_1 ) = J_EE_1( idx_J_EE_1 );
    J_EI( idx_J_EI_1 ) = J_EI_1( idx_J_EI_1 );

    for h1 = 1 : length( N_E )
        for h2 = 1 : length( N_E )
            if vars.connects( h1, h2 ) > 0
                J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) = J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) - idx_J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) .* ( sum( J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) - init_J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ), 2 ) ./ sum( idx_J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ), 2 ) );
            end
        end
    end
    J_EE( isnan( J_EE ) ) = 0;
    J_EE( idx_J_EE & J_EE < J_EE_min ) = J_EE_min;
    J_EE( idx_J_EE & J_EE > J_EE_max ) = J_EE_max;
end

% -------------------------------------------------------------------------

vars.J_EE = J_EE;
vars.J_EI = J_EI;

vars.V_E = V_E;
vars.V_I = V_I;
vars.Vu_E = Vu_E;
vars.Vv_E = Vv_E;

vars.s_E = s_E;
vars.sx_E = sx_E;
vars.sy_E = sy_E;
vars.s_I = s_I;
vars.sy_I = sy_I;

vars.V_T_E = V_T_E;
vars.ww_E = ww_E;

vars.ct_abs_E = ct_abs_E;
vars.ct_abs_I = ct_abs_I;

% vars.triplet_r_1 = triplet_r_1;
% vars.triplet_r_2 = triplet_r_2;
% vars.triplet_o_1 = triplet_o_1;
% vars.triplet_o_2 = triplet_o_2;

% -------------------------------------------------------------------------

% figure
% subplot( 2, 2, 1 )
% imagesc( V_E, [ -70, 30 ] )
% colorbar
% subplot( 2, 2, 2 )
% imagesc( V_I, [ -70, 30 ] )
% colorbar
% subplot( 2, 2, 2 + 1 )
% imagesc( s_E, [ 0, 1 ] )
% colorbar
% subplot( 2, 2, 2 + 2 )
% imagesc( s_I, [ 0, 1 ] )
% colorbar
