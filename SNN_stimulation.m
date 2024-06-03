function vars = SNN_stimulation( vars, lambda )

r_ext_EE = vars.r_ext_EE / 1000;
r_ext_IE = vars.r_ext_IE / 1000;

N_unit_whole = vars.N_E( vars.inputs );

nStims = vars.infoStims( 5 );

vars.ext_stim_E = zeros( sum( vars.N_E, 2 ), vars.infoStims( 1 ) / vars.dt );
vars.ext_stim_I = zeros( sum( vars.N_I, 2 ), vars.infoStims( 1 ) / vars.dt );
ct_E_h = 0;
for h = 1 : length( vars.inputs )
    if h == 1
        ext_stim_E = zeros( N_unit_whole( h ), vars.infoStims( 1 ) / vars.dt );
        tItv = 1 : ( vars.infoStims( 1 ) / vars.dt ) / round( ( vars.infoStims( 1 ) / vars.dt ) / ( vars.infoStims( 2 ) / vars.dt ) ) : ( vars.infoStims( 1 ) / vars.dt );
        nStimN = round( N_unit_whole( h ) / vars.infoStims( 4 ) );
        nItv = 1 : nStimN : N_unit_whole( h );
        ct_n = 0;
        for t = 1 : nStims
            ct_n = ct_n + 1;
            sn = mod( ct_n - 1, vars.infoStims( 4 ) ) + 1;
            ext_stim_E( nItv( sn ) + ( [ 1 : nStimN ] - 1 ), tItv( t ) + ( [ 1 : ( vars.infoStims( 1 ) / vars.dt ) / round( ( vars.infoStims( 1 ) / vars.dt ) / ( vars.infoStims( 2 ) / vars.dt ) ) ] - 1 ) ) = 1;
        end
    else
        ext_stim_E = circshift( ext_stim_E, ( vars.infoStims( 3 ) / vars.dt ), 2 );
    end
    if ismember( h, vars.inputs )
        vars.ext_stim_E( ct_E_h + [ 1 : vars.N_E( h ) ], : ) = ext_stim_E;
    end
    ct_E_h = ct_E_h + vars.N_E( h );
end

ext_noise_E = pinknoise( sum( vars.N_E, 2 ), vars.infoStims( 1 ) / vars.dt );
ext_noise_E = abs( ext_noise_E );
ext_noise_E( randperm( size( ext_noise_E, 1 ), round( 0.9 * size( ext_noise_E, 1 ) ) ), : ) = 0;
ext_noise_E = ext_noise_E / max( ext_noise_E, [], [ 1, 2 ] );
ext_noise_I = pinknoise( sum( vars.N_I, 2 ), vars.infoStims( 1 ) / vars.dt );
ext_noise_I = abs( ext_noise_I );
ext_noise_I( randperm( size( ext_noise_I, 1 ), round( 0.9 * size( ext_noise_I, 1 ) ) ), : ) = 0;
ext_noise_I = ext_noise_I / max( ext_noise_I, [], [ 1, 2 ] );

vars.ext_stim_E = vars.ext_stim_E + lambda * ext_noise_E;
vars.ext_stim_I = vars.ext_stim_I + lambda * ext_noise_I;

vars.ext_stim_E = r_ext_EE * vars.ext_stim_E;
vars.ext_stim_I = r_ext_IE * vars.ext_stim_I;
