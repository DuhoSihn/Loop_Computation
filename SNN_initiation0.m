function vars = SNN_initiation0( vars )

N_E = vars.N_E;
N_I = vars.N_I;

p_C = vars.p_C;
% p_E = vars.p_E;

J_EE_min = vars.J_EE_min;
J_EE_max = vars.J_EE_max;
J_EE_O = vars.J_EE_O;
J_EE_O1 = vars.J_EE_O1;
J_EI_min = vars.J_EI_min;
J_EI_max = vars.J_EI_max;
J_EI_O = vars.J_EI_O;
J_IE_val = vars.J_IE_val;
J_II_val = vars.J_II_val;
% J_GII_val = vars.J_GII_val;

nChannels = vars.nChannels;

% -------------------------------------------------------------------------

idxMat_E = false( length( N_E ), sum( N_E, 2 ) );
idxMatCh_E = false( length( N_E ), sum( N_E, 2 ), nChannels );
idxMatChAll_E = false( length( N_E ) * nChannels, sum( N_E, 2 ) );
ct1 = 0;
ct3 = 0;
for h = 1 : length( N_E )
    idxMat_E( h, ct1 + [ 1 : N_E( h ) ] ) = true;
    areaPartition = round( linspace( 0, N_E( h ), nChannels + 1 ) );
    for ch = 1 : nChannels
        idxMatCh_E( h, ct1 + [ areaPartition( ch ) + 1 : areaPartition( ch + 1 ) ], ch ) = true;
        ct3 = ct3 + 1;
        idxMatChAll_E( ct3, : ) = idxMatCh_E( h, :, ch );
    end
    ct1 = ct1 + N_E( h );
end
idxMat_I = false( length( N_I ), sum( N_I, 2 ) );
idxMatCh_I = false( length( N_I ), sum( N_I, 2 ), nChannels );
idxMatChAll_I = false( length( N_I ) * nChannels, sum( N_I, 2 ) );
ct1 = 0;
ct3 = 0;
for h = 1 : length( N_I )
    idxMat_I( h, ct1 + [ 1 : N_I( h ) ] ) = true;
    areaPartition = round( linspace( 0, N_I( h ), nChannels + 1 ) );
    for ch = 1 : nChannels
        idxMatCh_I( h, ct1 + [ areaPartition( ch ) + 1 : areaPartition( ch + 1 ) ], ch ) = true;
        ct3 = ct3 + 1;
        idxMatChAll_I( ct3, : ) = idxMatCh_I( h, :, ch );
    end
    ct1 = ct1 + N_I( h );
end

% -------------------------------------------------------------------------
% initiation of weights

idx_J_EE = false( sum( N_E, 2 ), sum( N_E, 2 ) );
idx_J_EI = false( sum( N_E, 2 ), sum( N_I, 2 ) );
idx_J_IE = false( sum( N_I, 2 ), sum( N_E, 2 ) );
idx_J_II = false( sum( N_I, 2 ), sum( N_I, 2 ) );

idx_J_EE_1 = false( sum( N_E, 2 ), sum( N_E, 2 ) );
idx_J_EI_1 = false( sum( N_E, 2 ), sum( N_I, 2 ) );
idx_J_IE_1 = false( sum( N_I, 2 ), sum( N_E, 2 ) );
idx_J_II_1 = false( sum( N_I, 2 ), sum( N_I, 2 ) );

J_EE = zeros( sum( N_E, 2 ), sum( N_E, 2 ) );
J_EI = zeros( sum( N_E, 2 ), sum( N_I, 2 ) );
J_IE = zeros( sum( N_I, 2 ), sum( N_E, 2 ) );
J_II = zeros( sum( N_I, 2 ), sum( N_I, 2 ) );

for h1 = 1 : length( N_E )
    for h2 = 1 : length( N_E )
        if vars.connects( h1, h2 ) > 0
            if h1 ~= h2
                idx_J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) = true;
                idx_J_EE_1( idxMat_E( h1, : ), idxMat_E( h2, : ) ) = true;
                %                 for ch = 1 : nChannels
                %                     idx_J_EE( idxMatCh_E( h1, :, ch ), idxMatCh_E( h2, :, ch ) ) = true;
                %                     idx_J_EE_1( idxMatCh_E( h1, :, ch ), idxMatCh_E( h2, :, ch ) ) = true;
                %                 end
            elseif h1 == h2
                % idx_J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) = true;
                if nChannels > 1
                    for ch = 1 : nChannels
                        idx_J_EE( idxMatCh_E( h1, :, ch ), idxMatCh_E( h2, :, ch ) ) = true;
                        % idx_J_EE( idxMatCh_E( h1, :, ch ), idxMatCh_E( h2, :, ch ) ) = false;
                        idx_J_EI( idxMatCh_E( h1, :, ch ), idxMatCh_I( h2, :, ch ) ) = true;
                        idx_J_EI_1( idxMatCh_E( h1, :, ch ), idxMatCh_I( h2, :, ch ) ) = true;
                        idx_J_IE( idxMatCh_I( h1, :, ch ), idxMatCh_E( h2, :, ch ) ) = true;
                        idx_J_IE_1( idxMatCh_I( h1, :, ch ), idxMatCh_E( h2, :, ch ) ) = true;
                        idx_J_II( idxMatCh_I( h1, :, ch ), idxMatCh_I( h2, :, ch ) ) = true;
                        idx_J_II_1( idxMatCh_I( h1, :, ch ), idxMatCh_I( h2, :, ch ) ) = true;
                    end
                elseif nChannels == 1
                    idx_J_EI( idxMat_E( h1, : ), idxMat_I( h2, : ) ) = true;
                    idx_J_EI_1( idxMat_E( h1, : ), idxMat_I( h2, : ) ) = true;
                    idx_J_IE( idxMat_I( h1, : ), idxMat_E( h2, : ) ) = true;
                    idx_J_IE_1( idxMat_I( h1, : ), idxMat_E( h2, : ) ) = true;
                    idx_J_II( idxMat_I( h1, : ), idxMat_I( h2, : ) ) = true;
                    idx_J_II_1( idxMat_I( h1, : ), idxMat_I( h2, : ) ) = true;
                end
            end
            if h1 ~= h2
                J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) = J_EE_O;
                %                 for ch = 1 : nChannels
                %                     J_EE( idxMatCh_E( h1, :, ch ), idxMatCh_E( h2, :, ch ) ) = J_EE_O;
                %                 end
                for ch = 1 : nChannels
                    if h2 == 1
                        J_EE( idxMatCh_E( h1, :, ch ), idxMatCh_E( h2, :, ch ) ) = J_EE_O1;
                    elseif h2 == 2
                        if ch < nChannels
                            J_EE( idxMatCh_E( h1, :, ch + 1 ), idxMatCh_E( h2, :, ch ) ) = J_EE_O1;
                        elseif ch == nChannels
                            J_EE( idxMatCh_E( h1, :, 1 ), idxMatCh_E( h2, :, ch ) ) = J_EE_O1;
                        end
                    end
                end
            elseif h1 == h2
                % J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) = J_EE_O;
                if nChannels > 1
                    for ch = 1 : nChannels
                        J_EE( idxMatCh_E( h1, :, ch ), idxMatCh_E( h2, :, ch ) ) = J_EE_O;
                        % J_EE( idxMatCh_E( h1, :, ch ), idxMatCh_E( h2, :, ch ) ) = 0;
                        J_EI( idxMatCh_E( h1, :, ch ), idxMatCh_I( h2, :, ch ) ) = J_EI_O;
                        J_IE( idxMatCh_I( h1, :, ch ), idxMatCh_E( h2, :, ch ) ) = J_IE_val;
                        J_II( idxMatCh_I( h1, :, ch ), idxMatCh_I( h2, :, ch ) ) = J_II_val;
                    end
                elseif nChannels == 1
                    J_EI( idxMat_E( h1, : ), idxMat_I( h2, : ) ) = J_EI_O;
                    J_IE( idxMat_I( h1, : ), idxMat_E( h2, : ) ) = J_IE_val;
                    J_II( idxMat_I( h1, : ), idxMat_I( h2, : ) ) = J_II_val;
                end
            end
        end
    end
end

p_idx = rand( size( idx_J_EE ) ) > p_C;
idx_J_EE( p_idx ) = false;
J_EE( p_idx ) = 0;

p_idx = rand( size( idx_J_EI ) ) > p_C;
idx_J_EI( p_idx ) = false;
J_EI( p_idx ) = 0;

p_idx = rand( size( idx_J_IE ) ) > p_C;
idx_J_IE( p_idx ) = false;
J_IE( p_idx ) = 0;

p_idx = rand( size( idx_J_II ) ) > p_C;
idx_J_II( p_idx ) = false;
J_II( p_idx ) = 0;

% -------------------------------------------------------------------------

scaleWeight = vars.scaleWeight;

fill_EE = mean( idx_J_EE( : ), 1 );
fill_EI = mean( idx_J_EI( : ), 1 );
fill_IE = mean( idx_J_IE( : ), 1 );
fill_II = mean( idx_J_II( : ), 1 );

J_EE = ( scaleWeight / fill_EE ) * J_EE;
J_EI = ( scaleWeight / fill_EI ) * J_EI;
J_IE = ( scaleWeight / fill_IE ) * J_IE;
J_II = ( scaleWeight / fill_II ) * J_II;

J_EE_min = ( scaleWeight / fill_EE ) * J_EE_min;
J_EE_max = ( scaleWeight / fill_EE ) * J_EE_max;
J_EE_O = ( scaleWeight / fill_EE ) * J_EE_O;
J_EI_min = ( scaleWeight / fill_EI ) * J_EI_min;
J_EI_max = ( scaleWeight / fill_EI ) * J_EI_max;

% -------------------------------------------------------------------------

init_J_EE = J_EE_O * idx_J_EE;

% Synaptic scaling
for h1 = 1 : length( N_E )
    for h2 = 1 : length( N_E )
        if vars.connects( h1, h2 ) > 0
            % J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) = J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) - idx_J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) .* ( sum( J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) - init_J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ), 2 ) ./ sum( idx_J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ), 2 ) );
            J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) = vars.connects( h1, h2 ) * ( J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) - idx_J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) .* ( sum( J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ) - init_J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ), 2 ) ./ sum( idx_J_EE( idxMat_E( h1, : ), idxMat_E( h2, : ) ), 2 ) ) );
        end
    end
end
% J_EE = J_EE - idx_J_EE .* ( sum( J_EE - init_J_EE, 2 ) ./ sum( idx_J_EE, 2 ) );
J_EE( isnan( J_EE ) ) = 0;
J_EE( idx_J_EE & J_EE < J_EE_min ) = J_EE_min;
J_EE( idx_J_EE & J_EE > J_EE_max ) = J_EE_max;

vars.idx_J_EE = idx_J_EE;
vars.idx_J_EI = idx_J_EI;
vars.idx_J_IE = idx_J_IE;
vars.idx_J_II = idx_J_II;
% vars.idx_J_GII = idx_J_GII;

vars.idx_J_EE_1 = idx_J_EE_1;
vars.idx_J_EI_1 = idx_J_EI_1;
vars.idx_J_IE_1 = idx_J_IE_1;
vars.idx_J_II_1 = idx_J_II_1;

vars.J_EE = J_EE;
vars.J_EI = J_EI;
vars.J_IE = J_IE;
vars.J_II = J_II;
% vars.J_GII = J_GII;

vars.init_J_EE = init_J_EE;
vars.idxMat_E = idxMat_E;

vars.J_EE_min = J_EE_min;
vars.J_EE_max = J_EE_max;
vars.J_EE_O = J_EE_O;
vars.J_EI_min = J_EI_min;
vars.J_EI_max = J_EI_max;

% -------------------------------------------------------------------------
%
% figure
%
% subplot( 2, 5, 1 )
% imagesc( idx_J_EE )
% subplot( 2, 5, 2 )
% imagesc( idx_J_EI )
% subplot( 2, 5, 3 )
% imagesc( idx_J_IE )
% subplot( 2, 5, 4 )
% imagesc( idx_J_II )
% subplot( 2, 5, 5 )
% imagesc( idx_J_GII )
%
% subplot( 2, 5, 5 + 1 )
% imagesc( J_EE )
% subplot( 2, 5, 5 + 2 )
% imagesc( J_EI )
% subplot( 2, 5, 5 + 3 )
% imagesc( J_IE )
% subplot( 2, 5, 5 + 4 )
% imagesc( J_II )
% subplot( 2, 5, 5 + 5 )
% imagesc( J_GII )
