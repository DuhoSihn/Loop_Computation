% Main_Loop_Computation
clear all; close all; clc
addpath(genpath('C:\DATA\Util\matlab_codes'));
addpath(genpath('D:\cDATA\Util\matlab_codes'));



%% Derive, Activity model
% %
% % T_exp = 200;
% % 
% % N_exc = 80;
% % N_inh = 20;
% % 
% % D_epsp = 19 / 20;
% % D_ipsp = 9 / 10;
% % 
% % threshold = 0.2;
% % threshold_pr = 0.5;
% % period_refractory = 3;
% % N_stim = 10;
% % 
% % W_ee = 1.0 * 2 / N_exc;
% % W_ei = 1.3 * 2 / N_exc;
% % W_ie = 1.3 * 2 / N_inh;
% % W_ii = 1.0 * 2 / N_inh;
% % 
% % W = [ W_ee * rand( N_exc, N_exc ), W_ie * rand( N_exc, N_inh ); W_ei * rand( N_inh, N_exc ), W_ii * rand( N_inh, N_inh ) ];
% % W( logical( eye( N_exc + N_inh ) ) ) = 0;
% % W_ext = [ 0.25 * W_ee * rand( N_exc, N_exc ), zeros( N_exc, N_inh ); zeros( N_inh, N_exc ), zeros( N_inh, N_inh ) ];
% % 
% % % figure
% % % subplot( 2, 1, 1 ); imagesc( W ); colorbar; subplot( 2, 1, 2 ); imagesc( W_ext ); colorbar;
% % 
% % 
% % count_refractory = zeros( N_exc + N_inh, 1 );
% % X_act = zeros( N_exc + N_inh, 1 );
% % X_psp = zeros( N_exc + N_inh, N_exc + N_inh );
% % 
% % r_act = zeros( N_exc + N_inh, T_exp );
% % r_psp = zeros( N_exc + N_inh, N_exc + N_inh, T_exp );
% % 
% % for t = 1 : T_exp
% % 
% %     % count refractory period
% %     count_refractory( count_refractory > 0 ) = count_refractory( count_refractory > 0 ) + 1;
% %     count_refractory( count_refractory > period_refractory ) = 0;
% % 
% %     % add external PSP
% %     X_ext = zeros( N_exc + N_inh, 1 );
% %     X_ext( randperm( N_exc, N_stim ), 1 ) = 1;
% %     X_psp = X_psp + bsxfun( @times, W_ext, transpose( X_ext ) );
% % 
% %     % add PSP
% %     X_psp = X_psp + bsxfun( @times, W, transpose( X_act ) );
% %     % leaky PSP
% %     X_psp = X_psp .* [ D_epsp * ones( N_exc + N_inh, N_exc ), D_ipsp * ones( N_exc + N_inh, N_inh ) ];
% %     X_psp( X_psp < 0 ) = 0;
% % 
% %     % action potential
% %     X_act = sum( X_psp( :, 1 : N_exc ), 2 ) - sum( X_psp( :, N_exc + 1 : N_exc + N_inh ), 2 ) > threshold & count_refractory == 0;
% %     X_act = X_act & rand( N_exc + N_inh, 1 ) > threshold_pr;
% %     count_refractory( X_act, 1 ) = count_refractory( X_act, 1 ) + 1;
% % 
% %     % recording
% %     r_act( :, t ) = X_act;
% %     r_psp( :, :, t ) = X_psp;
% % 
% % end; clear t
% % 
% % 
% % save( 'Results_model.mat', 'T_exp', 'N_exc', 'N_inh', 'D_epsp', 'D_ipsp', 'threshold', 'threshold_pr', 'period_refractory', 'N_stim', 'W_ee', 'W_ei', 'W_ie', 'W_ii', 'W', 'W_ext', 'r_act', 'r_psp' );
% %
% %
% load( 'Results_model.mat' )
% 
% 
% % epsp = 1;
% % ipsp = 1;
% % for t = 2 : T_exp
% %     epsp( t, 1 ) = D_epsp * epsp( t - 1, 1 );
% %     ipsp( t, 1 ) = D_ipsp * ipsp( t - 1, 1 );
% % end; clear t
% % figure
% % hold on
% % plot( [ 0 : T_exp - 1 ], epsp, 'r' )
% % plot( [ 0 : T_exp - 1 ], ipsp, 'b' );
% % clear epsp ipsp
% 
% 
% lw1 = 2;
% lw2 = 1.5;
% 
% figure( 'position', [ 100, 100, 700, 300 ] )
% subplot( 4, 1, [ 1, 2 ] )
% hold on
% plot( -100, 100, '-r', 'linewidth', lw1 )
% plot( -100, 100, '-b', 'linewidth', lw1 )
% for n = 1 : N_exc
%     if ~isempty( find( r_act( n, : ) == 1 ) )
%         plot( [ find( r_act( n, : ) == 1 ); find( r_act( n, : ) == 1 ) ], -n + [ 0.2 ; 0.8 ], '-r', 'linewidth', lw1 )
%     end
% end; clear n
% for n = N_exc + 1 : N_exc + N_inh
%     if ~isempty( find( r_act( n, : ) == 1 ) )
%         plot( [ find( r_act( n, : ) == 1 ); find( r_act( n, : ) == 1 ) ], -n + [ 0.2 ; 0.8 ], '-b', 'linewidth', lw1 )
%     end
% end; clear n
% set( gca, 'xlim', [ -5, 30 ], 'ylim', [ -( N_exc + N_inh ), 0 ], 'ytick', [] )
% ylabel( 'Neuron' )
% legend( 'Excitatory', 'Inhibitory', 'location', 'northwest' )
% subplot( 4, 1, 3 )
% hold on
% plot( smooth( 100 * sum( r_act( 1 : N_exc, : ), 1 ) / N_exc, 3 ), '-r', 'linewidth', lw2 )
% plot( smooth( 100 * sum( r_act( N_exc + 1 : N_exc + N_inh, : ), 1 ) / N_inh, 3 ), '-b', 'linewidth', lw2 )
% set( gca, 'xlim', [ -5, 30 ] )
% ylabel( 'Spike (%)' )
% subplot( 4, 1, 4 )
% hold on
% plot( permute( mean( mean( r_psp( 1 : N_exc, 1 : N_exc, : ), 2 ), 1 ), [ 3, 1, 2 ] ), '-r', 'linewidth', lw2 )
% plot( permute( mean( mean( r_psp( N_exc + 1 : N_exc + N_inh, 1 : N_exc, : ), 2 ), 1 ), [ 3, 1, 2 ] ), '-b', 'linewidth', lw2 )
% set( gca, 'xlim', [ -5, 30 ] )
% ylabel( 'PSP (a.u.)' )
% xlabel( 'Time (ms)' )
% 

%% Conditions for the emergence of sequential memory (functions of delta)
% 
% param_gamma = 50;% ms
% param_epsilon = 50;% ms
% param_delta = [ 0 : 50 ];% ms
% 
% N = [ 2 : 20 ];
% 
% 
% range_delta = param_gamma ./ ( N - 1 );
% range_delta_opt = param_gamma ./ N;
% 
% likelihood_f = ( ( 1 - ( transpose( param_delta ) ./ param_gamma ) ) .^ ( N - 1 ) ) .* ( ( ( N - 1 ) .* transpose( param_delta ) ) ./ param_gamma );
% likelihood_f( transpose( param_delta ) > range_delta ) = nan;
% likelihood_f_opt = ( 1 - ( 1 ./ N ) ) .^ N;% -> exp(-1) = 0.3679 as N -> inf
% 
% 
% lw = 2;
% colors1 = [ 1, 0, 1 ];
% colors2 = 0.2 * colors1 + 0.8 * [ 1, 1, 1 ];
% colors3 = [ 0, 0, 0 ];
% 
% 
% figure( 'position', [ 100, 100, 330, 300 ] )
% hold on
% patch( [ N, fliplr( N ) ], [ ones( 1, length( N ) ), fliplr( range_delta ) ], colors2, 'edgecolor', colors2 )
% plot( N, range_delta_opt, '-', 'color', colors1, 'linewidth', lw )
% set( gca, 'xlim', [ 1.9, max( N ) + 0.1 ], 'ylim', [ 0, max( param_delta ) ] )
% ylabel( 'Initiation period (\delta, ms)' )
% xlabel( 'Size of loop structure (N)' )
% title( 'Active and Inhibition periods (\gamma, \epsilon) = (50, 50) ms' )
% legend( 'Necessary condition', 'Optimal', 'Location', 'northeast')
% 
% 
% figure( 'position', [ 100, 100, 390, 300 ] )
% imagesc( N, param_delta, likelihood_f, 'alphadata', ~isnan( likelihood_f ), [ 0, max( likelihood_f( : ) ) ] )
% colormap( gca, turbo )
% colorbar
% axis xy
% ylabel( 'Initiation period (\delta, ms)' )
% xlabel( 'Size of loop structure (N)' )
% title( 'Likelihood (f(\delta))' )
% 
% 
% figure( 'position', [ 100, 100, 330, 300 ] )
% hold on
% plot( N, likelihood_f_opt, '-', 'color', colors3, 'linewidth', lw )
% plot( N, exp( -1 ) * ones( 1, length( N ) ), '--k' )
% set( gca, 'xlim', [ 1.9, max( N ) + 0.1 ], 'ylim', [ 0.24, 0.37 ] )
% ylabel( 'Optimal likelihood (max\{f(\delta)\})' )
% xlabel( 'Size of loop structure (N)' )
% 
% 
% % -------------------------------------------------------------------------
% 
% 
% param_gamma = 30;% ms
% param_epsilon = 30;% ms
% param_delta = [ 0 : 50 ];% ms
% 
% N = [ 2 : 20 ];
% 
% 
% range_delta = param_gamma ./ ( N - 1 );
% range_delta_opt = param_gamma ./ N;
% 
% likelihood_f = ( ( 1 - ( transpose( param_delta ) ./ param_gamma ) ) .^ ( N - 1 ) ) .* ( ( ( N - 1 ) .* transpose( param_delta ) ) ./ param_gamma );
% likelihood_f( transpose( param_delta ) > range_delta ) = nan;
% likelihood_f_opt = ( 1 - ( 1 ./ N ) ) .^ N;% -> exp(-1) = 0.3679 as N -> inf
% 
% 
% lw = 2;
% colors1 = [ 1, 0, 1 ];
% colors2 = 0.2 * colors1 + 0.8 * [ 1, 1, 1 ];
% colors3 = [ 0, 0, 0 ];
% 
% 
% figure( 'position', [ 100, 100, 330, 300 ] )
% hold on
% patch( [ N, fliplr( N ) ], [ ones( 1, length( N ) ), fliplr( range_delta ) ], colors2, 'edgecolor', colors2 )
% plot( N, range_delta_opt, '-', 'color', colors1, 'linewidth', lw )
% set( gca, 'xlim', [ 1.9, max( N ) + 0.1 ], 'ylim', [ 0, max( param_delta ) ] )
% ylabel( 'Initiation period (\delta, ms)' )
% xlabel( 'Size of loop structure (N)' )
% title( 'Active and Inhibition periods (\gamma, \epsilon) = (30, 30) ms' )
% legend( 'Necessary condition', 'Optimal', 'Location', 'northeast')
% 
% 
% figure( 'position', [ 100, 100, 390, 300 ] )
% imagesc( N, param_delta, likelihood_f, 'alphadata', ~isnan( likelihood_f ), [ 0, max( likelihood_f( : ) ) ] )
% colormap( gca, turbo )
% colorbar
% axis xy
% ylabel( 'Initiation period (\delta, ms)' )
% xlabel( 'Size of loop structure (N)' )
% title( 'Likelihood (f(\delta))' )
% 
% 
% figure( 'position', [ 100, 100, 330, 300 ] )
% hold on
% plot( N, likelihood_f_opt, '-', 'color', colors3, 'linewidth', lw )
% plot( N, exp( -1 ) * ones( 1, length( N ) ), '--k' )
% set( gca, 'xlim', [ 1.9, max( N ) + 0.1 ], 'ylim', [ 0.24, 0.37 ] )
% ylabel( 'Optimal likelihood (max\{f(\delta)\})' )
% xlabel( 'Size of loop structure (N)' )
% 

%% Conditions for the emergence of sequential memory (functions of gamma)
% 
% param_gamma = [ 0 : 50 ];% ms
% param_epsilon = 50;% ms
% param_delta = 5;% ms
% 
% N = [ 2 : 20 ];
% 
% 
% range_gamma = param_delta .* ( N - 1 );
% range_gamma_opt = param_delta .* N;
% 
% likelihood_f = ( ( 1 - ( param_delta ./ transpose( param_gamma ) ) ) .^ ( N - 1 ) ) .* ( ( ( N - 1 ) .* param_delta ) ./ transpose( param_gamma ) );
% likelihood_f( transpose( param_gamma ) < range_gamma ) = nan;
% likelihood_f_opt = ( 1 - ( 1 ./ N ) ) .^ N;% -> exp(-1) = 0.3679 as N -> inf
% 
% 
% lw = 2;
% colors1 = [ 0, 1, 1 ];
% colors2 = 0.2 * colors1 + 0.8 * [ 1, 1, 1 ];
% colors3 = [ 0, 0, 0 ];
% 
% 
% figure( 'position', [ 100, 100, 330, 300 ] )
% hold on
% patch( [ N, fliplr( N ) ], [ 50.5 * ones( 1, length( N ) ), fliplr( range_gamma ) ], colors2, 'edgecolor', colors2 )
% plot( N, range_gamma_opt, '-', 'color', colors1, 'linewidth', lw )
% set( gca, 'xlim', [ 1.9, max( N ) + 0.1 ], 'ylim', [ 0, max( param_gamma ) ] )
% ylabel( 'Active period (\gamma, ms)' )
% xlabel( 'Size of loop structure (N)' )
% title( 'Initiation periods (\delta) = 5 ms' )
% legend( 'Necessary condition', 'Optimal', 'Location', 'southeast')
% 
% 
% figure( 'position', [ 100, 100, 390, 300 ] )
% imagesc( N, param_gamma, likelihood_f, 'alphadata', ~isnan( likelihood_f ), [ 0, max( likelihood_f( : ) ) ] )
% colormap( gca, turbo )
% colorbar
% axis xy
% ylabel( 'Active period (\gamma, ms)' )
% xlabel( 'Size of loop structure (N)' )
% title( 'Likelihood (f(\gamma))' )
% 
% 
% figure( 'position', [ 100, 100, 330, 300 ] )
% hold on
% plot( N, likelihood_f_opt, '-', 'color', colors3, 'linewidth', lw )
% plot( N, exp( -1 ) * ones( 1, length( N ) ), '--k' )
% set( gca, 'xlim', [ 1.9, max( N ) + 0.1 ], 'ylim', [ 0.24, 0.37 ] )
% ylabel( 'Optimal likelihood (max\{f(\gamma)\})' )
% xlabel( 'Size of loop structure (N)' )
% 

%% Get W1 and W2 for Example test
%
% pool_period_initiation = [ 0 : 30 ];
%
% Ws = cell( length( pool_period_initiation ), 2 );
%
% pop_size = 5;
% n_stim = 10;
% n_iter_learning = 1000;
%
% ct_period_initiation = 0;
% for period_initiation = pool_period_initiation
%     ct_period_initiation = ct_period_initiation + 1;
%
%     period_active = 60 - period_initiation;
%
%     S1 = [];
%     for k = 1 : n_stim
%         tS = [ zeros( n_stim * pop_size, period_active ) ];
%         tS( ( k - 1 ) * pop_size + [ 1 : pop_size ], : ) = 1;
%         S1 = [ S1, tS ];
%     end; clear k tS
%     S1 = repmat( S1, [ 1, n_iter_learning ] );
%     S2 = circshift( S1, period_initiation, 2);
%
%     tic
%     [ W1, W2 ] = learningHebbian( S1, S2 );
%     toc
%
%     Ws{ ct_period_initiation, 1 } = W1;
%     Ws{ ct_period_initiation, 2 } = W2;
%
%     clear period_active S1 S2 W1 W2
%
%     disp( [ num2str( period_initiation ), ' / ', num2str( pool_period_initiation( end ) ), ' parameters.' ] )
%
% end; clear period_initiation ct_period_initiation
%
%
% save( 'W1_W2_example.mat', 'pool_period_initiation', 'Ws', 'pop_size', 'n_stim', 'n_iter_learning' )
%

%% Example test
% 
% load( 'W1_W2_example.mat' )
% 
% threshold = 0.15;
% 
% S = zeros( pop_size * n_stim, n_stim );
% for k = 1 : n_stim
%     S( ( k - 1 ) * pop_size + [ 1 : pop_size ], k ) = 1;
% end; clear k
% 
% pool_period_refractory = [ 0 : 60 ];
% 
% M = [];
% C = [];
% H = [];
% 
% ct_period_initiation = 0;
% for period_initiation = pool_period_initiation
%     ct_period_initiation = ct_period_initiation + 1;
%     
%     period_active = 60 - period_initiation;
%     
%     ct_period_refractory = 0;
%     for period_refractory = pool_period_refractory
%         ct_period_refractory = ct_period_refractory + 1;
%         
%         S1 = [];
%         for k = 1 : n_stim
%             tS = [ zeros( n_stim * pop_size, period_active ) ];
%             tS( ( k - 1 ) * pop_size + [ 1 : pop_size ], : ) = 1;
%             S1 = [ S1, tS ];
%         end; clear k tS
%         S2 = circshift( S1, period_initiation, 2);
%         
%         W1 = Ws{ ct_period_initiation, 1 };
%         W2 = Ws{ ct_period_initiation, 2 };
%         
%         ct_tau = 0;
%         for tau = [ 0 : period_active : size( S1, 2 ) - 1 ]
%             ct_tau = ct_tau + 1;
%             
%             tS1 = [ circshift( S1, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS2 = [ circshift( S2, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             
%             [ Z1, Z2 ] = assemblyComputation( tS1, tS2, W1, W2, threshold, period_initiation, period_active, period_refractory );
%             
%             M1 = mean( sum( Z1( :, size( S1, 2 ) + 1 :  end ), 1 ), 2 );
%             M2 = mean( sum( Z2( :, size( S2, 2 ) + 1 :  end ), 1 ), 2 );
%             [ ~, C1, H1 ] = patternMatching( S, Z1( :, size( S1, 2 ) + 1 :  end ) );
%             [ ~, C2, H2 ] = patternMatching( S, Z2( :, size( S2, 2 ) + 1 :  end ) );
%             C1 = C1( H1 > 0 );
%             C2 = C2( H2 > 0 );
%             if length( C1 ) >= 0.5 * length( H1 )
%                 t1C = C1( 1 : end - 1 );
%                 t2C = C1( 2 : end );
%                 idxC = [];
%                 for n = 1 : n_stim
%                     if n < n_stim
%                         idx = t1C == n;
%                         idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                     elseif n == n_stim
%                         idx = t1C == n;
%                         idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                     end
%                     if ~isempty( idx )
%                         idxC = [ idxC, all( idx, 2 ) ];
%                     else
%                         idxC = [ idxC, false ];
%                     end
%                     clear idx
%                 end; clear n
%                 C1 = mean( idxC, 2 ) >= 0.5;
%                 clear t1C t2C idxC
%             else
%                 C1 = false;
%             end
%             if length( C2 ) >= 0.5 * length( H2 )
%                 t1C = C2( 1 : end - 1 );
%                 t2C = C2( 2 : end );
%                 idxC = [];
%                 for n = 1 : n_stim
%                     if n < n_stim
%                         idx = t1C == n;
%                         idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                     elseif n == n_stim
%                         idx = t1C == n;
%                         idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                     end
%                     if ~isempty( idx )
%                         idxC = [ idxC, all( idx, 2 ) ];
%                     else
%                         idxC = [ idxC, false ];
%                     end
%                     clear idx
%                 end; clear n
%                 C2 = mean( idxC, 2 ) >= 0.5;
%                 clear t1C t2C idxC
%             else
%                 C2 = false;
%             end
%             H1 = mean( H1, 2 );
%             H2 = mean( H2, 2 );
%             
%             M( ct_period_initiation, ct_period_refractory, ct_tau ) = mean( [ M1, M2 ], 2 );
%             C( ct_period_initiation, ct_period_refractory, ct_tau ) = all( [ C1, C2 ], 2 );
%             H( ct_period_initiation, ct_period_refractory, ct_tau ) = mean( [ H1, H2 ], 2 );
%             
%             clear tS1 tS2 Z1 Z2 M1 M2 C1 C2 H1 H2
%         end; clear tau ct_tau
%         
%         clear S1 S2 W1 W2
%         
%         disp( [ num2str( period_initiation ), ' / ', num2str( pool_period_initiation( end ) ), ' and ', num2str( period_refractory ), ' / ', num2str( pool_period_refractory( end ) ), ' parameters.' ] )
%         
%     end; clear period_refractory
%     
%     clear period_active
%     
% end; clear period_initiation
% 
% 
% save( 'Results_example.mat', 'M', 'C', 'H', 'threshold', 'pool_period_refractory', 'pool_period_initiation', 'Ws', 'pop_size', 'n_stim', 'n_iter_learning' )
% 

%% Results for Example test
% 
% load( 'Results_example.mat' )
% 
% 
% % -------------------------------------------------------------------------
% 
% 
% likelihood_f = ( 1 - ( transpose( pool_period_initiation ) ./ ( 60 - transpose( pool_period_initiation ) ) ) ) .* ( transpose( pool_period_initiation ) ./ ( 60 - transpose( pool_period_initiation ) ) );
% likelihood_f = repmat( likelihood_f, [ 1, length( pool_period_refractory ) ] );
% 
% likelihood_f( repmat( transpose( pool_period_initiation ), [ 1, length( pool_period_refractory ) ] ) <= 0 ) = NaN;
% likelihood_f( bsxfun( @gt, transpose( pool_period_initiation ), 60 - pool_period_initiation ) ) = NaN;
% likelihood_f( bsxfun( @gt, transpose( pool_period_initiation ), pool_period_refractory ) ) = NaN;
% 
% 
% figure( 'position', [ 100, 100, 300, 200 ] )
% imagesc( pool_period_refractory, pool_period_initiation, likelihood_f, 'alphadata', ~isnan( likelihood_f ), [ 0, 0.25 ] )
% colormap( gca, turbo )
% colorbar
% axis xy
% ylabel( 'Initiation period (ms)' )
% xlabel( 'Inhibition period (ms)' )
% title( 'Theoretical likelihood' )
% 
% 
% % -------------------------------------------------------------------------
% 
% 
% colors = [ 1, 0.3, 0; 0, 1, 0.3; 0.3, 0, 1; 0.8, 0.8, 0.3 ];
% 
% 
% M = mean( M, 3 );
% C = all( C, 3 );
% H = mean( H, 3 );
% 
% P = NaN( size( C ) );
% P( M < 1 ) = 1;
% P( M > 15 ) = 2;
% P( M >= 1 & M <= 15 & C == 0 ) = 3;
% P( M >= 1 & M <= 15 & C == 1 ) = 4;
% 
% % figure( 'position', [ 100, 100, 300, 200 ] )
% % imagesc( pool_period_refractory, pool_period_initiation, M, [ 0, max( M( : ), [], 1 ) ] )
% % colormap( gca, turbo )
% % colorbar
% % ylabel( 'Initiation period (ms)' )
% % xlabel( 'Inhibition period (ms)' )
% % title( 'Number of activated assemblies' )
% 
% % figure( 'position', [ 100, 100, 300, 200 ] )
% % imagesc( pool_period_refractory, pool_period_initiation, C, [ 0, 1 ] )
% % colormap( gca, [ 1, 1, 1; colors( 4, : ) ] )
% % colorbar
% % ylabel( 'Initiation period (ms)' )
% % xlabel( 'Inhibition period (ms)' )
% % title( 'Sequence generation' )
% 
% figure( 'position', [ 100, 100, 300, 200 ] )
% imagesc( pool_period_refractory, pool_period_initiation, P, [ 1, 4 ] )
% colormap( gca, colors )
% colorbar
% axis xy
% ylabel( 'Initiation period (ms)' )
% xlabel( 'Inhibition period (ms)' )
% title( 'States of parameters' )
% 
% 
% % -------------------------------------------------------------------------
% 
% 
% ct_period_initiation = 1;
% ct_period_refractory = 41;
% period_initiation = pool_period_initiation( ct_period_initiation );
% period_active = 60 - period_initiation;
% period_refractory = pool_period_refractory( ct_period_refractory );
% 
% S1 = [];
% for k = 1 : n_stim
%     tS = [ zeros( n_stim * pop_size, period_active ) ];
%     tS( ( k - 1 ) * pop_size + [ 1 : pop_size ], : ) = 1;
%     S1 = [ S1, tS ];
% end; clear k tS
% S2 = circshift( S1, period_initiation, 2);
% 
% W1 = Ws{ ct_period_initiation, 1 };
% W2 = Ws{ ct_period_initiation, 2 };
% 
% tau = 0;
% tS1 = [ circshift( S1, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% tS2 = [ circshift( S2, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% 
% [ Z1, Z2 ] = assemblyComputation( tS1, tS2, W1, W2, threshold, period_initiation, period_active, period_refractory );
% 
% figure( 'position', [ 100, 100, 500, 120 ] )
% imagesc( Z1 )
% colormap( gca, [ 1, 1, 1; colors( 1, : ) ] )
% hold on
% plot( ( size( S1, 2 ) + 0.5 ) * [ 1, 1 ], [ 0.5, pop_size * n_stim + 0.5 ], '-k' )
% set( gca, 'xlim', [ 0.5, 3000.5 ], 'ylim', [ 0.5, pop_size * n_stim + 0.5 ], 'ytick', [] )
% ylabel( 'Assembly' )
% xlabel( 'Time (ms)' )
% title( [ 'Passive (Initiation period: ', num2str( period_initiation ), ' ms / Inhibition period: ', num2str( period_refractory ), ' ms )' ] )
% 
% clear ct_period_initiation ct_period_refractory period_initiation period_active period_refractory S1 S2 W1 W2 tau tS1 tS2 Z1 Z2
% 
% 
% % -------------------------------------------------------------------------
% 
% 
% ct_period_initiation = 21;
% ct_period_refractory = 1;
% period_initiation = pool_period_initiation( ct_period_initiation );
% period_active = 60 - period_initiation;
% period_refractory = pool_period_refractory( ct_period_refractory );
% 
% S1 = [];
% for k = 1 : n_stim
%     tS = [ zeros( n_stim * pop_size, period_active ) ];
%     tS( ( k - 1 ) * pop_size + [ 1 : pop_size ], : ) = 1;
%     S1 = [ S1, tS ];
% end; clear k tS
% S2 = circshift( S1, period_initiation, 2);
% 
% W1 = Ws{ ct_period_initiation, 1 };
% W2 = Ws{ ct_period_initiation, 2 };
% 
% tau = 0;
% tS1 = [ circshift( S1, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% tS2 = [ circshift( S2, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% 
% [ Z1, Z2 ] = assemblyComputation( tS1, tS2, W1, W2, threshold, period_initiation, period_active, period_refractory );
% 
% figure( 'position', [ 100, 100, 500, 120 ] )
% imagesc( Z1 )
% colormap( gca, [ 1, 1, 1; colors( 2, : ) ] )
% hold on
% plot( ( size( S1, 2 ) + 0.5 ) * [ 1, 1 ], [ 0.5, pop_size * n_stim + 0.5 ], '-k' )
% set( gca, 'xlim', [ 0.5, 3000.5 ], 'ylim', [ 0.5, pop_size * n_stim + 0.5 ], 'ytick', [] )
% ylabel( 'Assembly' )
% xlabel( 'Time (ms)' )
% title( [ 'Seizure-like (Initiation period: ', num2str( period_initiation ), ' ms / Inhibition period: ', num2str( period_refractory ), ' ms )' ] )
% 
% clear ct_period_initiation ct_period_refractory period_initiation period_active period_refractory S1 S2 W1 W2 tau tS1 tS2 Z1 Z2
% 
% 
% % -------------------------------------------------------------------------
% 
% 
% ct_period_initiation = 1;
% ct_period_refractory = 1;
% period_initiation = pool_period_initiation( ct_period_initiation );
% period_active = 60 - period_initiation;
% period_refractory = pool_period_refractory( ct_period_refractory );
% 
% S1 = [];
% for k = 1 : n_stim
%     tS = [ zeros( n_stim * pop_size, period_active ) ];
%     tS( ( k - 1 ) * pop_size + [ 1 : pop_size ], : ) = 1;
%     S1 = [ S1, tS ];
% end; clear k tS
% S2 = circshift( S1, period_initiation, 2);
% 
% W1 = Ws{ ct_period_initiation, 1 };
% W2 = Ws{ ct_period_initiation, 2 };
% 
% tau = 0;
% tS1 = [ circshift( S1, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% tS2 = [ circshift( S2, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% 
% [ Z1, Z2 ] = assemblyComputation( tS1, tS2, W1, W2, threshold, period_initiation, period_active, period_refractory );
% 
% figure( 'position', [ 100, 100, 500, 120 ] )
% imagesc( Z1 )
% colormap( gca, [ 1, 1, 1; colors( 3, : ) ] )
% hold on
% plot( ( size( S1, 2 ) + 0.5 ) * [ 1, 1 ], [ 0.5, pop_size * n_stim + 0.5 ], '-k' )
% set( gca, 'xlim', [ 0.5, 3000.5 ], 'ylim', [ 0.5, pop_size * n_stim + 0.5 ], 'ytick', [] )
% ylabel( 'Assembly' )
% xlabel( 'Time (ms)' )
% title( [ 'Persistent (Initiation period: ', num2str( period_initiation ), ' ms / Inhibition period: ', num2str( period_refractory ), ' ms )' ] )
% 
% clear ct_period_initiation ct_period_refractory period_initiation period_active period_refractory S1 S2 W1 W2 tau tS1 tS2 Z1 Z2
% 
% 
% % -------------------------------------------------------------------------
% 
% 
% ct_period_initiation = 21;
% ct_period_refractory = 41;
% period_initiation = pool_period_initiation( ct_period_initiation );
% period_active = 60 - period_initiation;
% period_refractory = pool_period_refractory( ct_period_refractory );
% 
% S1 = [];
% for k = 1 : n_stim
%     tS = [ zeros( n_stim * pop_size, period_active ) ];
%     tS( ( k - 1 ) * pop_size + [ 1 : pop_size ], : ) = 1;
%     S1 = [ S1, tS ];
% end; clear k tS
% S2 = circshift( S1, period_initiation, 2);
% 
% W1 = Ws{ ct_period_initiation, 1 };
% W2 = Ws{ ct_period_initiation, 2 };
% 
% tau = 0;
% tS1 = [ circshift( S1, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% tS2 = [ circshift( S2, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% 
% [ Z1, Z2 ] = assemblyComputation( tS1, tS2, W1, W2, threshold, period_initiation, period_active, period_refractory );
% 
% figure( 'position', [ 100, 100, 500, 120 ] )
% imagesc( Z1 )
% colormap( gca, [ 1, 1, 1; colors( 4, : ) ] )
% hold on
% plot( ( size( S1, 2 ) + 0.5 ) * [ 1, 1 ], [ 0.5, pop_size * n_stim + 0.5 ], '-k' )
% set( gca, 'xlim', [ 0.5, 3000.5 ], 'ylim', [ 0.5, pop_size * n_stim + 0.5 ], 'ytick', [] )
% ylabel( 'Assembly' )
% xlabel( 'Time (ms)' )
% title( [ 'Self-generation (Initiation period: ', num2str( period_initiation ), ' ms / Inhibition period: ', num2str( period_refractory ), ' ms )' ] )
% 
% clear ct_period_initiation ct_period_refractory period_initiation period_active period_refractory S1 S2 W1 W2 tau tS1 tS2 Z1 Z2
% 
% 
% % -------------------------------------------------------------------------
% 
% 
% ct_period_initiation = 21;
% ct_period_refractory = 41;
% period_initiation = pool_period_initiation( ct_period_initiation );
% period_active = 60 - period_initiation;
% period_refractory = pool_period_refractory( ct_period_refractory );
% 
% S1 = [];
% for k = 1 : n_stim
%     tS = [ zeros( n_stim * pop_size, period_active ) ];
%     tS( ( k - 1 ) * pop_size + [ 1 : pop_size ], : ) = 1;
%     S1 = [ S1, tS ];
% end; clear k tS
% S2 = circshift( S1, period_initiation, 2);
% 
% 
% figure( 'position', [ 100, 100, 300, 150 ] )
% subplot( 1, 2, 1 )
% imagesc( S1 )
% colormap( gca, [ 1, 1, 1; 0, 0, 0 ] )
% set( gca, 'ytick', [] )
% ylabel( 'Assembly' )
% xlabel( 'Time (ms)' )
% title( 'Area X' )
% subplot( 1, 2, 2 )
% imagesc( S2 )
% colormap( gca, [ 1, 1, 1; 0, 0, 0 ] )
% set( gca, 'ytick', [] )
% xlabel( 'Time (ms)' )
% title( 'Area Y' )
% 

%% Get W1 and W2 for Example test (N = 3)
% 
% pool_period_initiation = [ 0 : 30 ];
% 
% Ws = cell( length( pool_period_initiation ), 3 );
% 
% pop_size = 5;
% n_stim = 10;
% n_iter_learning = 1000;
% 
% ct_period_initiation = 0;
% for period_initiation = pool_period_initiation
%     ct_period_initiation = ct_period_initiation + 1;
% 
%     period_active = 60 - period_initiation;
% 
%     S1 = [];
%     for k = 1 : n_stim
%         tS = [ zeros( n_stim * pop_size, period_active ) ];
%         tS( ( k - 1 ) * pop_size + [ 1 : pop_size ], : ) = 1;
%         S1 = [ S1, tS ];
%     end; clear k tS
%     S1 = repmat( S1, [ 1, n_iter_learning ] );
%     S2 = circshift( S1, period_initiation, 2);
%     S3 = circshift( S2, period_initiation, 2);
% 
%     tic
%     [ W1, W2, W3 ] = learningHebbian3( S1, S2, S3 );
%     toc
% 
%     Ws{ ct_period_initiation, 1 } = W1;
%     Ws{ ct_period_initiation, 2 } = W2;
%     Ws{ ct_period_initiation, 3 } = W3;
% 
%     clear period_active S1 S2 S3 W1 W2 W3
% 
%     disp( [ num2str( period_initiation ), ' / ', num2str( pool_period_initiation( end ) ), ' parameters.' ] )
% 
% end; clear period_initiation ct_period_initiation
% 
% 
% save( 'W1_W2_W3_example.mat', 'pool_period_initiation', 'Ws', 'pop_size', 'n_stim', 'n_iter_learning' )
% 

%% Example test (N = 3)
% 
% load( 'W1_W2_W3_example.mat' )
% 
% threshold = 0.15;
% 
% S = zeros( pop_size * n_stim, n_stim );
% for k = 1 : n_stim
%     S( ( k - 1 ) * pop_size + [ 1 : pop_size ], k ) = 1;
% end; clear k
% 
% pool_period_refractory = [ 0 : 60 ];
% 
% M = [];
% C = [];
% H = [];
% 
% ct_period_initiation = 0;
% for period_initiation = pool_period_initiation
%     ct_period_initiation = ct_period_initiation + 1;
%     
%     period_active = 60 - period_initiation;
%     
%     ct_period_refractory = 0;
%     for period_refractory = pool_period_refractory
%         ct_period_refractory = ct_period_refractory + 1;
%         
%         S1 = [];
%         for k = 1 : n_stim
%             tS = [ zeros( n_stim * pop_size, period_active ) ];
%             tS( ( k - 1 ) * pop_size + [ 1 : pop_size ], : ) = 1;
%             S1 = [ S1, tS ];
%         end; clear k tS
%         S2 = circshift( S1, period_initiation, 2);
%         S3 = circshift( S2, period_initiation, 2);
%         
%         W1 = Ws{ ct_period_initiation, 1 };
%         W2 = Ws{ ct_period_initiation, 2 };
%         W3 = Ws{ ct_period_initiation, 3 };
%         
%         ct_tau = 0;
%         for tau = [ 0 : period_active : size( S1, 2 ) - 1 ]
%             ct_tau = ct_tau + 1;
%             
%             tS1 = [ circshift( S1, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS2 = [ circshift( S2, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS3 = [ circshift( S3, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             
%             [ Z1, Z2, Z3 ] = assemblyComputation3( tS1, tS2, tS3, W1, W2, W3, threshold, period_initiation, period_active, period_refractory );
%             
%             M1 = mean( sum( Z1( :, size( S1, 2 ) + 1 :  end ), 1 ), 2 );
%             M2 = mean( sum( Z2( :, size( S2, 2 ) + 1 :  end ), 1 ), 2 );
%             M3 = mean( sum( Z3( :, size( S3, 2 ) + 1 :  end ), 1 ), 2 );
%             [ ~, C1, H1 ] = patternMatching( S, Z1( :, size( S1, 2 ) + 1 :  end ) );
%             [ ~, C2, H2 ] = patternMatching( S, Z2( :, size( S2, 2 ) + 1 :  end ) );
%             [ ~, C3, H3 ] = patternMatching( S, Z3( :, size( S3, 2 ) + 1 :  end ) );
%             C1 = C1( H1 > 0 );
%             C2 = C2( H2 > 0 );
%             C3 = C3( H3 > 0 );
%             if length( C1 ) >= 0.5 * length( H1 )
%                 t1C = C1( 1 : end - 1 );
%                 t2C = C1( 2 : end );
%                 idxC = [];
%                 for n = 1 : n_stim
%                     if n < n_stim
%                         idx = t1C == n;
%                         idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                     elseif n == n_stim
%                         idx = t1C == n;
%                         idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                     end
%                     if ~isempty( idx )
%                         idxC = [ idxC, all( idx, 2 ) ];
%                     else
%                         idxC = [ idxC, false ];
%                     end
%                     clear idx
%                 end; clear n
%                 C1 = mean( idxC, 2 ) >= 0.5;
%                 clear t1C t2C idxC
%             else
%                 C1 = false;
%             end
%             if length( C2 ) >= 0.5 * length( H2 )
%                 t1C = C2( 1 : end - 1 );
%                 t2C = C2( 2 : end );
%                 idxC = [];
%                 for n = 1 : n_stim
%                     if n < n_stim
%                         idx = t1C == n;
%                         idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                     elseif n == n_stim
%                         idx = t1C == n;
%                         idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                     end
%                     if ~isempty( idx )
%                         idxC = [ idxC, all( idx, 2 ) ];
%                     else
%                         idxC = [ idxC, false ];
%                     end
%                     clear idx
%                 end; clear n
%                 C2 = mean( idxC, 2 ) >= 0.5;
%                 clear t1C t2C idxC
%             else
%                 C2 = false;
%             end
%             if length( C3 ) >= 0.5 * length( H3 )
%                 t1C = C3( 1 : end - 1 );
%                 t2C = C3( 2 : end );
%                 idxC = [];
%                 for n = 1 : n_stim
%                     if n < n_stim
%                         idx = t1C == n;
%                         idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                     elseif n == n_stim
%                         idx = t1C == n;
%                         idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                     end
%                     if ~isempty( idx )
%                         idxC = [ idxC, all( idx, 2 ) ];
%                     else
%                         idxC = [ idxC, false ];
%                     end
%                     clear idx
%                 end; clear n
%                 C3 = mean( idxC, 2 ) >= 0.5;
%                 clear t1C t2C idxC
%             else
%                 C3 = false;
%             end
%             H1 = mean( H1, 2 );
%             H2 = mean( H2, 2 );
%             H3 = mean( H3, 2 );
%             
%             M( ct_period_initiation, ct_period_refractory, ct_tau ) = mean( [ M1, M2, M3 ], 2 );
%             C( ct_period_initiation, ct_period_refractory, ct_tau ) = all( [ C1, C2, C2 ], 2 );
%             H( ct_period_initiation, ct_period_refractory, ct_tau ) = mean( [ H1, H2, H3 ], 2 );
%             
%             clear tS1 tS2 tS3 Z1 Z2 Z3 M1 M2 M3 C1 C2 C3 H1 H2 H3
%         end; clear tau ct_tau
%         
%         clear S1 S2 S3 W1 W2 W3
%         
%         disp( [ num2str( period_initiation ), ' / ', num2str( pool_period_initiation( end ) ), ' and ', num2str( period_refractory ), ' / ', num2str( pool_period_refractory( end ) ), ' parameters.' ] )
%         
%     end; clear period_refractory
%     
%     clear period_active
%     
% end; clear period_initiation
% 
% 
% save( 'Results_example3.mat', 'M', 'C', 'H', 'threshold', 'pool_period_refractory', 'pool_period_initiation', 'Ws', 'pop_size', 'n_stim', 'n_iter_learning' )
% 

%% Results for Example test (N = 3)
% 
% load( 'Results_example3.mat' )
% 
% 
% % -------------------------------------------------------------------------
% 
% N = 3;
% 
% likelihood_f = ( ( 1 - ( transpose( pool_period_initiation ) ./ ( 60 - transpose( pool_period_initiation ) ) ) ) .^ ( N - 1 ) ) .* ( ( N - 1 ) * transpose( pool_period_initiation ) ./ ( 60 - transpose( pool_period_initiation ) ) );
% likelihood_f = repmat( likelihood_f, [ 1, length( pool_period_refractory ) ] );
% 
% likelihood_f( repmat( transpose( pool_period_initiation ), [ 1, length( pool_period_refractory ) ] ) <= 0 ) = NaN;
% likelihood_f( bsxfun( @gt, transpose( pool_period_initiation ), ( 60 - pool_period_initiation ) / ( N - 1 ) ) ) = NaN;
% likelihood_f( bsxfun( @gt, transpose( pool_period_initiation ), pool_period_refractory / ( N - 1 ) ) ) = NaN;
% 
% 
% figure( 'position', [ 100, 100, 300, 200 ] )
% imagesc( pool_period_refractory, pool_period_initiation, likelihood_f, 'alphadata', ~isnan( likelihood_f ), [ 0, 0.30 ] )
% colormap( gca, turbo )
% colorbar
% axis xy
% ylabel( 'Initiation period (ms)' )
% xlabel( 'Inhibition period (ms)' )
% title( 'Theoretical likelihood' )
% 
% 
% % -------------------------------------------------------------------------
% 
% 
% colors = [ 1, 0.3, 0; 0, 1, 0.3; 0.3, 0, 1; 0.8, 0.8, 0.3 ];
% 
% 
% M = mean( M, 3 );
% C = all( C, 3 );
% H = mean( H, 3 );
% 
% P = NaN( size( C ) );
% P( M < 1 ) = 1;
% P( M > 15 ) = 2;
% P( M >= 1 & M <= 15 & C == 0 ) = 3;
% P( M >= 1 & M <= 15 & C == 1 ) = 4;
% 
% % figure( 'position', [ 100, 100, 300, 200 ] )
% % imagesc( pool_period_refractory, pool_period_initiation, M, [ 0, max( M( : ), [], 1 ) ] )
% % colormap( gca, turbo )
% % colorbar
% % axis xy
% % ylabel( 'Initiation period (ms)' )
% % xlabel( 'Inhibition period (ms)' )
% % title( 'Number of activated assemblies' )
% 
% % figure( 'position', [ 100, 100, 300, 200 ] )
% % imagesc( pool_period_refractory, pool_period_initiation, C, [ 0, 1 ] )
% % colormap( gca, [ 1, 1, 1; colors( 4, : ) ] )
% % colorbar
% % axis xy
% % ylabel( 'Initiation period (ms)' )
% % xlabel( 'Inhibition period (ms)' )
% % title( 'Sequence generation' )
% 
% figure( 'position', [ 100, 100, 300, 200 ] )
% imagesc( pool_period_refractory, pool_period_initiation, P, [ 1, 4 ] )
% colormap( gca, colors )
% colorbar
% axis xy
% hold on
% plot( 30, 0, 'ok' )
% plot( 10, 10, 'ok' )
% plot( 35, 10, 'ok' )
% plot( 0, 2, 'ok' )
% plot( 35, 16, 'ok' )
% ylabel( 'Initiation period (ms)' )
% xlabel( 'Inhibition period (ms)' )
% title( 'States of parameters' )
% 
% 
% % -------------------------------------------------------------------------
% 
% 
% ct_period_initiation = 0 + 1;
% ct_period_refractory = 30 + 1;
% period_initiation = pool_period_initiation( ct_period_initiation );
% period_active = 60 - period_initiation;
% period_refractory = pool_period_refractory( ct_period_refractory );
% 
% S1 = [];
% for k = 1 : n_stim
%     tS = [ zeros( n_stim * pop_size, period_active ) ];
%     tS( ( k - 1 ) * pop_size + [ 1 : pop_size ], : ) = 1;
%     S1 = [ S1, tS ];
% end; clear k tS
% S2 = circshift( S1, period_initiation, 2);
% S3 = circshift( S2, period_initiation, 2);
% 
% W1 = Ws{ ct_period_initiation, 1 };
% W2 = Ws{ ct_period_initiation, 2 };
% W3 = Ws{ ct_period_initiation, 3 };
% 
% tau = 0;
% tS1 = [ circshift( S1, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% tS2 = [ circshift( S2, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% tS3 = [ circshift( S3, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% 
% [ Z1, Z2, Z3 ] = assemblyComputation3( tS1, tS2, tS3, W1, W2, W3, threshold, period_initiation, period_active, period_refractory );
% 
% figure( 'position', [ 100, 100, 500, 120 ] )
% imagesc( Z1 )
% colormap( gca, [ 1, 1, 1; colors( 1, : ) ] )
% hold on
% plot( ( size( S1, 2 ) + 0.5 ) * [ 1, 1 ], [ 0.5, pop_size * n_stim + 0.5 ], '-k' )
% set( gca, 'xlim', [ 0.5, 3000.5 ], 'ylim', [ 0.5, pop_size * n_stim + 0.5 ], 'ytick', [] )
% ylabel( 'Assembly' )
% xlabel( ' ' )
% % xlabel( 'Time (ms)' )
% title( [ 'Passive (Initiation period: ', num2str( period_initiation ), ' ms / Inhibition period: ', num2str( period_refractory ), ' ms )' ] )
% 
% clear ct_period_initiation ct_period_refractory period_initiation period_active period_refractory S1 S2 W1 W2 tau tS1 tS2 Z1 Z2
% 
% 
% % -------------------------------------------------------------------------
% 
% 
% ct_period_initiation = 10 + 1;
% ct_period_refractory = 10 + 1;
% period_initiation = pool_period_initiation( ct_period_initiation );
% period_active = 60 - period_initiation;
% period_refractory = pool_period_refractory( ct_period_refractory );
% 
% S1 = [];
% for k = 1 : n_stim
%     tS = [ zeros( n_stim * pop_size, period_active ) ];
%     tS( ( k - 1 ) * pop_size + [ 1 : pop_size ], : ) = 1;
%     S1 = [ S1, tS ];
% end; clear k tS
% S2 = circshift( S1, period_initiation, 2);
% S3 = circshift( S2, period_initiation, 2);
% 
% W1 = Ws{ ct_period_initiation, 1 };
% W2 = Ws{ ct_period_initiation, 2 };
% W3 = Ws{ ct_period_initiation, 3 };
% 
% tau = 0;
% tS1 = [ circshift( S1, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% tS2 = [ circshift( S2, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% tS3 = [ circshift( S3, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% 
% [ Z1, Z2, Z3 ] = assemblyComputation3( tS1, tS2, tS3, W1, W2, W3, threshold, period_initiation, period_active, period_refractory );
% 
% figure( 'position', [ 100, 100, 500, 120 ] )
% imagesc( Z1 )
% colormap( gca, [ 1, 1, 1; colors( 2, : ) ] )
% hold on
% plot( ( size( S1, 2 ) + 0.5 ) * [ 1, 1 ], [ 0.5, pop_size * n_stim + 0.5 ], '-k' )
% set( gca, 'xlim', [ 0.5, 3000.5 ], 'ylim', [ 0.5, pop_size * n_stim + 0.5 ], 'ytick', [] )
% ylabel( 'Assembly' )
% xlabel( ' ' )
% % xlabel( 'Time (ms)' )
% title( [ 'Seizure-like (Initiation period: ', num2str( period_initiation ), ' ms / Inhibition period: ', num2str( period_refractory ), ' ms )' ] )
% 
% clear ct_period_initiation ct_period_refractory period_initiation period_active period_refractory S1 S2 W1 W2 tau tS1 tS2 Z1 Z2
% 
% 
% % -------------------------------------------------------------------------
% 
% 
% ct_period_initiation = 10 + 1;
% ct_period_refractory = 35 + 1;
% period_initiation = pool_period_initiation( ct_period_initiation );
% period_active = 60 - period_initiation;
% period_refractory = pool_period_refractory( ct_period_refractory );
% 
% S1 = [];
% for k = 1 : n_stim
%     tS = [ zeros( n_stim * pop_size, period_active ) ];
%     tS( ( k - 1 ) * pop_size + [ 1 : pop_size ], : ) = 1;
%     S1 = [ S1, tS ];
% end; clear k tS
% S2 = circshift( S1, period_initiation, 2);
% S3 = circshift( S2, period_initiation, 2);
% 
% W1 = Ws{ ct_period_initiation, 1 };
% W2 = Ws{ ct_period_initiation, 2 };
% W3 = Ws{ ct_period_initiation, 3 };
% 
% tau = 0;
% tS1 = [ circshift( S1, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% tS2 = [ circshift( S2, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% tS3 = [ circshift( S3, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% 
% [ Z1, Z2, Z3 ] = assemblyComputation3( tS1, tS2, tS3, W1, W2, W3, threshold, period_initiation, period_active, period_refractory );
% 
% figure( 'position', [ 100, 100, 500, 120 ] )
% imagesc( Z1 )
% colormap( gca, [ 1, 1, 1; colors( 2, : ) ] )
% hold on
% plot( ( size( S1, 2 ) + 0.5 ) * [ 1, 1 ], [ 0.5, pop_size * n_stim + 0.5 ], '-k' )
% set( gca, 'xlim', [ 0.5, 3000.5 ], 'ylim', [ 0.5, pop_size * n_stim + 0.5 ], 'ytick', [] )
% ylabel( 'Assembly' )
% xlabel( ' ' )
% % xlabel( 'Time (ms)' )
% title( [ 'Seizure-like (Initiation period: ', num2str( period_initiation ), ' ms / Inhibition period: ', num2str( period_refractory ), ' ms )' ] )
% 
% clear ct_period_initiation ct_period_refractory period_initiation period_active period_refractory S1 S2 W1 W2 tau tS1 tS2 Z1 Z2
% 
% 
% % -------------------------------------------------------------------------
% 
% 
% ct_period_initiation = 2 + 1;
% ct_period_refractory = 0 + 1;
% period_initiation = pool_period_initiation( ct_period_initiation );
% period_active = 60 - period_initiation;
% period_refractory = pool_period_refractory( ct_period_refractory );
% 
% S1 = [];
% for k = 1 : n_stim
%     tS = [ zeros( n_stim * pop_size, period_active ) ];
%     tS( ( k - 1 ) * pop_size + [ 1 : pop_size ], : ) = 1;
%     S1 = [ S1, tS ];
% end; clear k tS
% S2 = circshift( S1, period_initiation, 2);
% S3 = circshift( S2, period_initiation, 2);
% 
% W1 = Ws{ ct_period_initiation, 1 };
% W2 = Ws{ ct_period_initiation, 2 };
% W3 = Ws{ ct_period_initiation, 3 };
% 
% tau = 0;
% tS1 = [ circshift( S1, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% tS2 = [ circshift( S2, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% tS3 = [ circshift( S3, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% 
% [ Z1, Z2, Z3 ] = assemblyComputation3( tS1, tS2, tS3, W1, W2, W3, threshold, period_initiation, period_active, period_refractory );
% 
% figure( 'position', [ 100, 100, 500, 120 ] )
% imagesc( Z1 )
% colormap( gca, [ 1, 1, 1; colors( 3, : ) ] )
% hold on
% plot( ( size( S1, 2 ) + 0.5 ) * [ 1, 1 ], [ 0.5, pop_size * n_stim + 0.5 ], '-k' )
% set( gca, 'xlim', [ 0.5, 3000.5 ], 'ylim', [ 0.5, pop_size * n_stim + 0.5 ], 'ytick', [] )
% ylabel( 'Assembly' )
% xlabel( ' ' )
% % xlabel( 'Time (ms)' )
% title( [ 'Persistent (Initiation period: ', num2str( period_initiation ), ' ms / Inhibition period: ', num2str( period_refractory ), ' ms )' ] )
% 
% clear ct_period_initiation ct_period_refractory period_initiation period_active period_refractory S1 S2 W1 W2 tau tS1 tS2 Z1 Z2
% 
% 
% % -------------------------------------------------------------------------
% 
% 
% ct_period_initiation = 16 + 1;
% ct_period_refractory = 35 + 1;
% period_initiation = pool_period_initiation( ct_period_initiation );
% period_active = 60 - period_initiation;
% period_refractory = pool_period_refractory( ct_period_refractory );
% 
% S1 = [];
% for k = 1 : n_stim
%     tS = [ zeros( n_stim * pop_size, period_active ) ];
%     tS( ( k - 1 ) * pop_size + [ 1 : pop_size ], : ) = 1;
%     S1 = [ S1, tS ];
% end; clear k tS
% S2 = circshift( S1, period_initiation, 2);
% S3 = circshift( S2, period_initiation, 2);
% 
% W1 = Ws{ ct_period_initiation, 1 };
% W2 = Ws{ ct_period_initiation, 2 };
% W3 = Ws{ ct_period_initiation, 3 };
% 
% tau = 0;
% tS1 = [ circshift( S1, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% tS2 = [ circshift( S2, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% tS3 = [ circshift( S3, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
% 
% [ Z1, Z2, Z3 ] = assemblyComputation3( tS1, tS2, tS3, W1, W2, W3, threshold, period_initiation, period_active, period_refractory );
% 
% figure( 'position', [ 100, 100, 500, 120 ] )
% imagesc( Z1 )
% colormap( gca, [ 1, 1, 1; colors( 4, : ) ] )
% hold on
% plot( ( size( S1, 2 ) + 0.5 ) * [ 1, 1 ], [ 0.5, pop_size * n_stim + 0.5 ], '-k' )
% set( gca, 'xlim', [ 0.5, 3000.5 ], 'ylim', [ 0.5, pop_size * n_stim + 0.5 ], 'ytick', [] )
% ylabel( 'Assembly' )
% xlabel( 'Time (ms)' )
% title( [ 'Self-generation (Initiation period: ', num2str( period_initiation ), ' ms / Inhibition period: ', num2str( period_refractory ), ' ms )' ] )
% 
% clear ct_period_initiation ct_period_refractory period_initiation period_active period_refractory S1 S2 W1 W2 tau tS1 tS2 Z1 Z2
% 
% 
% % % -------------------------------------------------------------------------
% 
% 
% ct_period_initiation = 16 + 1;
% ct_period_refractory = 35 + 1;
% period_initiation = pool_period_initiation( ct_period_initiation );
% period_active = 60 - period_initiation;
% period_refractory = pool_period_refractory( ct_period_refractory );
% 
% S1 = [];
% for k = 1 : n_stim
%     tS = [ zeros( n_stim * pop_size, period_active ) ];
%     tS( ( k - 1 ) * pop_size + [ 1 : pop_size ], : ) = 1;
%     S1 = [ S1, tS ];
% end; clear k tS
% S2 = circshift( S1, period_initiation, 2);
% S3 = circshift( S2, period_initiation, 2);
% 
% 
% figure( 'position', [ 100, 100, 300, 150 ] )
% subplot( 1, 3, 1 )
% imagesc( S1 )
% colormap( gca, [ 1, 1, 1; 0, 0, 0 ] )
% set( gca, 'ytick', [] )
% ylabel( 'Assembly' )
% xlabel( 'Time (ms)' )
% title( 'Area X' )
% subplot( 1, 3, 2 )
% imagesc( S2 )
% colormap( gca, [ 1, 1, 1; 0, 0, 0 ] )
% set( gca, 'ytick', [] )
% xlabel( 'Time (ms)' )
% title( 'Area Y' )
% subplot( 1, 3, 3 )
% imagesc( S3 )
% colormap( gca, [ 1, 1, 1; 0, 0, 0 ] )
% set( gca, 'ytick', [] )
% xlabel( 'Time (ms)' )
% title( 'Area Z' )
% 
