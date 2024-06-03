% Main_Loop_Computation
clear all; close all; clc



%% Derive, Basic neural activity model
%
% T_exp = 200;
%
% N_exc = 2 * 80;
% N_inh = 2 * 20;
%
% % dt = 1;
% dt = 0.1;
% t_memory = 100;
% t_memory_L = t_memory / dt;
% tau_c = 20;% membran potential time constant (ms)
% tau_r_E = 1 * 1;% rise time for E synapses (ms)
% tau_d_E = 6 * 1;% decay time for E synapses (ms)
% tau_r_I = 0.5 * 1;% rise time for I synapses (ms)
% tau_d_I = 2 * 1;% decay time for I synapses (ms)
% t_domain_memory = 0 : dt : t_memory - dt;
% F_E = ( 1 / ( tau_d_E - tau_r_E ) ) * ( exp( -t_domain_memory / tau_d_E ) - exp( -t_domain_memory / tau_r_E ) );
% F_I = ( 1 / ( tau_d_I - tau_r_I ) ) * ( exp( -t_domain_memory / tau_d_I ) - exp( -t_domain_memory / tau_r_I ) );
% % figure; hold on; plot(F_E); plot(F_I); set( gca, 'xlim', [ 0, length( t_domain_memory ) ] )
% F = fliplr( [ repmat( F_E, [ N_exc, 1 ] ); repmat( F_I, [ N_inh, 1 ] ) ] );
%
% % threshold = 0.2;
% threshold = 0.13;
% threshold_pr = 1 - 0.05 * dt;
% period_refractory = 3;
% N_stim = 20;
%
% W_ee = 0.4 * 2 / N_exc;
% W_ei = 1.6 * 2 / N_exc;
% W_ie = 1.6 * 2 / N_inh;
% W_ii = 1.6 * 2 / N_inh;
% % W_ee = 1.0 * 2 / N_exc;
% % W_ei = 1.3 * 2 / N_exc;
% % W_ie = 1.3 * 2 / N_inh;
% % W_ii = 1.0 * 2 / N_inh;
% % W_ee = 0.7 * 2 / N_exc;
% % W_ei = 1.3 * 2 / N_exc;
% % W_ie = 1.3 * 2 / N_inh;
% % W_ii = 0.7 * 2 / N_inh;
% % W_ee = 1.0 * 2 / N_exc;
% % W_ei = 1.0 * 2 / N_exc;
% % W_ie = 1.0 * 2 / N_inh;
% % W_ii = 1.0 * 2 / N_inh;
%
% thr_conn = 0.2;
% J_EE = rand( N_exc, N_exc );
% J_EI = rand( N_exc, N_inh );
% J_IE = rand( N_inh, N_exc );
% J_II = rand( N_inh, N_inh );
% J_EE_ext = rand( N_exc, N_exc );
% J_EE( rand( N_exc, N_exc ) > thr_conn ) = 0;
% J_EI( rand( N_exc, N_inh ) > thr_conn ) = 0;
% J_IE( rand( N_inh, N_exc ) > thr_conn ) = 0;
% J_II( rand( N_inh, N_inh ) > thr_conn ) = 0;
% J_EE_ext( rand( N_exc, N_exc ) > thr_conn ) = 0;
% % W = [ W_ee * rand( N_exc, N_exc ), W_ie * rand( N_exc, N_inh ); W_ei * rand( N_inh, N_exc ), W_ii * rand( N_inh, N_inh ) ];
% W = [ W_ee * J_EE, W_ie * J_EI; W_ei * J_IE, W_ii * J_II ];
% W( logical( eye( N_exc + N_inh ) ) ) = 0;
% % W_ext = [ 0.25 * W_ee * rand( N_exc, N_exc ), zeros( N_exc, N_inh ); zeros( N_inh, N_exc ), zeros( N_inh, N_inh ) ];
% % W_ext = [ 0.25 * W_ee * J_EE_ext, zeros( N_exc, N_inh ); zeros( N_inh, N_exc ), zeros( N_inh, N_inh ) ];
% W_ext = [ 0.7 * W_ee * J_EE_ext, zeros( N_exc, N_inh ); zeros( N_inh, N_exc ), zeros( N_inh, N_inh ) ];
%
% % figure
% % subplot( 2, 1, 1 ); imagesc( W ); colorbar; subplot( 2, 1, 2 ); imagesc( W_ext ); colorbar;
%
%
% count_refractory = zeros( N_exc + N_inh, 1 );
% X_act = zeros( N_exc + N_inh, 1 );
% X_psp = zeros( N_exc + N_inh, N_exc + N_inh );
%
% t_act = zeros( N_exc + N_inh, t_memory / dt );
%
% r_act = zeros( N_exc + N_inh, T_exp / dt );
% r_psp = zeros( N_exc + N_inh, N_exc + N_inh, T_exp / dt );
%
% ct_t = 0;
% for t = 0 : dt : T_exp
%     ct_t = ct_t + 1;
%
%     % count refractory period
%     count_refractory( count_refractory > 0 ) = count_refractory( count_refractory > 0 ) + 1 * dt;
%     count_refractory( count_refractory > period_refractory ) = 0;
%
%     % leaky PSP
%     X_psp = X_psp - X_psp / ( tau_c / dt );
%
%     % add external PSP
%     X_ext = zeros( N_exc + N_inh, 1 );
%     X_ext( randperm( N_exc, N_stim ), 1 ) = 1;
%     X_psp = X_psp + dt * bsxfun( @times, W_ext, transpose( X_ext ) );
%
%     % add PSP
%     % X_psp = X_psp + bsxfun( @times, W, transpose( X_act ) );
%     X_psp = X_psp + bsxfun( @times, W, transpose( sum( ( t_act .* F ), 2 ) ) );
%     X_psp( X_psp < 0 ) = 0;
%
%     % action potential
%     X_act = sum( X_psp( :, 1 : N_exc ), 2 ) - sum( X_psp( :, N_exc + 1 : N_exc + N_inh ), 2 ) > threshold & count_refractory == 0;
%     X_act = X_act & rand( N_exc + N_inh, 1 ) > threshold_pr;
%     count_refractory( X_act, 1 ) = count_refractory( X_act, 1 ) + 1 * dt;
%
%     % recording
%     r_act( :, ct_t ) = X_act;
%     r_psp( :, :, ct_t ) = X_psp;
%
%     % temporal recording
%     t_act = [ t_act( :, 2 : end ), X_act ];
%
% end; clear t ct_t
%
%
% save( 'Results_model.mat', 'T_exp', 'N_exc', 'N_inh', 'F', 'threshold', 'threshold_pr', 'period_refractory', 'N_stim', 'W_ee', 'W_ei', 'W_ie', 'W_ii', 'W', 'W_ext', 'r_act', 'r_psp' );
%
%
% load( 'Results_model.mat' )
%
%
% lw1 = 2;
% lw2 = 1.5;
% xlim2 = 80 / dt;
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
% set( gca, 'xlim', [ -5 / dt, xlim2 ], 'ylim', [ -( N_exc + N_inh ), 0 ], 'ytick', [] )
% set( gca, 'xtick', [ 0 : 10 : 80 ] / dt, 'xticklabel', { '0', '10', '20', '30', '40', '50', '60', '70', '80' } )
% ylabel( 'Neuron' )
% legend( 'Excitatory', 'Inhibitory', 'location', 'northwest' )
% subplot( 4, 1, 3 )
% hold on
% plot( smooth( 1000 * sum( r_act( 1 : N_exc, : ), 1 ) / N_exc, 3 / dt ), '-r', 'linewidth', lw2 )
% plot( smooth( 1000 * sum( r_act( N_exc + 1 : N_exc + N_inh, : ), 1 ) / N_inh, 3 / dt ), '-b', 'linewidth', lw2 )
% set( gca, 'xlim', [ -5 / dt, xlim2 ] )
% set( gca, 'xtick', [ 0 : 10 : 80 ] / dt, 'xticklabel', { '0', '10', '20', '30', '40', '50', '60', '70', '80' } )
% ylabel( 'Spikes/s' )
% subplot( 4, 1, 4 )
% hold on
% % plot( permute( mean( mean( r_psp( 1 : N_exc, 1 : N_exc, : ), 2 ), 1 ), [ 3, 1, 2 ] ), '-r', 'linewidth', lw2 )
% % plot( permute( mean( mean( r_psp( N_exc + 1 : N_exc + N_inh, 1 : N_exc, : ), 2 ), 1 ), [ 3, 1, 2 ] ), '-b', 'linewidth', lw2 )
% plot( permute( mean( sum( r_psp( 1 : N_exc, 1 : N_exc, : ), 2 ) - sum( r_psp( 1 : N_exc, N_exc + 1 : N_exc + N_inh, : ), 2 ), 1 ), [ 3, 1, 2 ] ), '-r', 'linewidth', lw2 )
% plot( permute( mean( sum( r_psp( N_exc + 1 : N_exc + N_inh, 1 : N_exc, : ), 2 ) - sum( r_psp( N_exc + 1 : N_exc + N_inh, N_exc + 1 : N_exc + N_inh, : ), 2 ), 1 ), [ 3, 1, 2 ] ), '-b', 'linewidth', lw2 )
% set( gca, 'xlim', [ -5 / dt, xlim2 ] )
% set( gca, 'xtick', [ 0 : 10 : 80 ] / dt, 'xticklabel', { '0', '10', '20', '30', '40', '50', '60', '70', '80' } )
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
% figure( 'position', [ 100, 100, 220, 200 ] )
% hold on
% patch( [ N, fliplr( N ) ], [ ones( 1, length( N ) ), fliplr( range_delta ) ], colors2, 'edgecolor', colors2 )
% plot( N, range_delta_opt, '-', 'color', colors1, 'linewidth', lw )
% set( gca, 'xlim', [ 1.9, max( N ) + 0.1 ], 'ylim', [ 0, max( param_delta ) ] )
% % ylabel( 'Initiation period (\delta, ms)' )
% % xlabel( 'Size of loop structure (N)' )
% % title( 'Active and Inhibition periods (\gamma, \epsilon) = (50, 50) ms' )
% legend( 'Necessary condition', 'Optimal', 'Location', 'northeast')
%
%
% figure( 'position', [ 100, 100, 260, 200 ] )
% imagesc( N, param_delta, likelihood_f, 'alphadata', ~isnan( likelihood_f ), [ 0, max( likelihood_f( : ) ) ] )
% colormap( gca, turbo )
% colorbar
% axis xy
% % ylabel( 'Initiation period (\delta, ms)' )
% % xlabel( 'Size of loop structure (N)' )
% % title( 'Likelihood (f(\delta))' )
%
%
% figure( 'position', [ 100, 100, 220, 200 ] )
% hold on
% plot( N, likelihood_f_opt, '-', 'color', colors3, 'linewidth', lw )
% plot( N, exp( -1 ) * ones( 1, length( N ) ), '--k' )
% set( gca, 'xlim', [ 1.9, max( N ) + 0.1 ], 'ylim', [ 0.24, 0.37 ] )
% % ylabel( 'Optimal likelihood (max\{f(\delta)\})' )
% % xlabel( 'Size of loop structure (N)' )
%
%
% % -------------------------------------------------------------------------
%
%
% param_gamma = 100;% ms
% param_epsilon = 100;% ms
% param_delta = [ 0 : 100 ];% ms
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
% figure( 'position', [ 100, 100, 220, 200 ] )
% hold on
% patch( [ N, fliplr( N ) ], [ ones( 1, length( N ) ), fliplr( range_delta ) ], colors2, 'edgecolor', colors2 )
% plot( N, range_delta_opt, '-', 'color', colors1, 'linewidth', lw )
% set( gca, 'xlim', [ 1.9, max( N ) + 0.1 ], 'ylim', [ 0, max( param_delta ) ] )
% % ylabel( 'Initiation period (\delta, ms)' )
% % xlabel( 'Size of loop structure (N)' )
% % title( 'Active and Inhibition periods (\gamma, \epsilon) = (80, 80) ms' )
% legend( 'Necessary condition', 'Optimal', 'Location', 'northeast')
%
%
% figure( 'position', [ 100, 100, 260, 200 ] )
% imagesc( N, param_delta, likelihood_f, 'alphadata', ~isnan( likelihood_f ), [ 0, max( likelihood_f( : ) ) ] )
% colormap( gca, turbo )
% colorbar
% axis xy
% % ylabel( 'Initiation period (\delta, ms)' )
% % xlabel( 'Size of loop structure (N)' )
% % title( 'Likelihood (f(\delta))' )
%
%
% figure( 'position', [ 100, 100, 220, 200 ] )
% hold on
% plot( N, likelihood_f_opt, '-', 'color', colors3, 'linewidth', lw )
% plot( N, exp( -1 ) * ones( 1, length( N ) ), '--k' )
% set( gca, 'xlim', [ 1.9, max( N ) + 0.1 ], 'ylim', [ 0.24, 0.37 ] )
% % ylabel( 'Optimal likelihood (max\{f(\delta)\})' )
% % xlabel( 'Size of loop structure (N)' )
%

%% Conditions for the emergence of sequential memory (functions of gamma)
%
% param_gamma = [ 0 : 50 ];% ms
% param_epsilon = 20;% ms
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
% figure( 'position', [ 100, 100, 220, 200 ] )
% hold on
% patch( [ N, fliplr( N ) ], [ 50.5 * ones( 1, length( N ) ), fliplr( range_gamma ) ], colors2, 'edgecolor', colors2 )
% plot( N, range_gamma_opt, '-', 'color', colors1, 'linewidth', lw )
% set( gca, 'xlim', [ 1.9, max( N ) + 0.1 ], 'ylim', [ 0, max( param_gamma ) ] )
% % ylabel( 'Active period (\gamma, ms)' )
% % xlabel( 'Size of loop structure (N)' )
% % title( 'Initiation periods (\delta) = 5 ms' )
% legend( 'Necessary condition', 'Optimal', 'Location', 'southeast')
%
%
% figure( 'position', [ 100, 100, 260, 200 ] )
% imagesc( N, param_gamma, likelihood_f, 'alphadata', ~isnan( likelihood_f ), [ 0, max( likelihood_f( : ) ) ] )
% colormap( gca, turbo )
% colorbar
% axis xy
% % ylabel( 'Active period (\gamma, ms)' )
% % xlabel( 'Size of loop structure (N)' )
% % title( 'Likelihood (f(\gamma))' )
%
%
% figure( 'position', [ 100, 100, 220, 200 ] )
% hold on
% plot( N, likelihood_f_opt, '-', 'color', colors3, 'linewidth', lw )
% plot( N, exp( -1 ) * ones( 1, length( N ) ), '--k' )
% set( gca, 'xlim', [ 1.9, max( N ) + 0.1 ], 'ylim', [ 0.24, 0.37 ] )
% % ylabel( 'Optimal likelihood (max\{f(\gamma)\})' )
% % xlabel( 'Size of loop structure (N)' )
%

%% Conditions for the emergence of sequential memory (functions of delta and gamma)
%
% param_gamma = [ 30 : 90 ];% ms
% param_epsilon = 50;% ms
% param_delta = [ 0 : 60 ];% ms
%
% N = 3;
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
% figure( 'position', [ 100, 100, 260, 200 ] )
% imagesc( param_delta, param_gamma, likelihood_f, 'alphadata', ~isnan( likelihood_f ), [ 0, max( likelihood_f( : ) ) ] )
% colormap( gca, turbo )
% colorbar
% axis xy
% % ylabel( 'Active period (\gamma, ms)' )
% % xlabel( 'Size of loop structure (N)' )
% % title( 'Likelihood (f(\gamma))' )
%

%% Get W1 and W2 for Example test (N = 2)
%
% pool_period_initiation = [ 20 : 60 ];
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
%     period_active = 120 - period_initiation;
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
%     [ W1, W2 ] = learningHebbian2( S1, S2 );
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

%% Example test (N = 2)
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
% pool_period_refractory = [ 0 : 120 ];
%
% M = [];
% C = [];
% H = [];
%
% ct_period_initiation = 0;
% for period_initiation = pool_period_initiation
%     ct_period_initiation = ct_period_initiation + 1;
%
%     period_active = 120 - period_initiation;
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
%             [ Z1, Z2 ] = assemblyComputation2( tS1, tS2, W1, W2, threshold, period_initiation, period_active, period_refractory );
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
% save( 'Results_example2.mat', 'M', 'C', 'H', 'threshold', 'pool_period_refractory', 'pool_period_initiation', 'Ws', 'pop_size', 'n_stim', 'n_iter_learning' )
%

%% Results for Example test (N = 2)
%
% load( 'Results_example2.mat' )
%
%
% % -------------------------------------------------------------------------
%
%
% likelihood_f = ( 1 - ( transpose( pool_period_initiation ) ./ ( 120 - transpose( pool_period_initiation ) ) ) ) .* ( transpose( pool_period_initiation ) ./ ( 120 - transpose( pool_period_initiation ) ) );
% likelihood_f = repmat( likelihood_f, [ 1, length( pool_period_refractory ) ] );
%
% likelihood_f( repmat( transpose( pool_period_initiation ), [ 1, length( pool_period_refractory ) ] ) <= 0 ) = NaN;
% likelihood_f( bsxfun( @gt, transpose( pool_period_initiation ), 120 - pool_period_initiation ) ) = NaN;
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

%% Get W1 and W2 for Example test (N = 3)
%
% pool_period_initiation = [ 20 : 60 ];
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
%     period_active = 120 - period_initiation;
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
% pool_period_refractory = [ 0 : 120 ];
%
% M = [];
% C = [];
% H = [];
%
% ct_period_initiation = 0;
% for period_initiation = pool_period_initiation
%     ct_period_initiation = ct_period_initiation + 1;
%
%     period_active = 120 - period_initiation;
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
%             C( ct_period_initiation, ct_period_refractory, ct_tau ) = all( [ C1, C2, C3 ], 2 );
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
% likelihood_f = ( ( 1 - ( transpose( pool_period_initiation ) ./ ( 120 - transpose( pool_period_initiation ) ) ) ) .^ ( N - 1 ) ) .* ( ( N - 1 ) * transpose( pool_period_initiation ) ./ ( 120 - transpose( pool_period_initiation ) ) );
% likelihood_f = repmat( likelihood_f, [ 1, length( pool_period_refractory ) ] );
%
% likelihood_f( repmat( transpose( pool_period_initiation ), [ 1, length( pool_period_refractory ) ] ) <= 0 ) = NaN;
% likelihood_f( bsxfun( @gt, transpose( pool_period_initiation ), ( 120 - pool_period_initiation ) / ( N - 1 ) ) ) = NaN;
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
% plot( 100, 60, 'ok' )
% plot( 20, 40, 'ok' )
% plot( 20, 58, 'ok' )
% plot( 100, 35, 'ok' )
% ylabel( 'Initiation period (ms)' )
% xlabel( 'Inhibition period (ms)' )
% title( 'States of parameters' )
%
%
% % -------------------------------------------------------------------------
%
%
% ct_period_initiation = 40 + 1;
% ct_period_refractory = 100 + 1;
% period_initiation = pool_period_initiation( ct_period_initiation );
% period_active = 120 - period_initiation;
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
% ct_period_initiation = 20 + 1;
% ct_period_refractory = 20 + 1;
% period_initiation = pool_period_initiation( ct_period_initiation );
% period_active = 120 - period_initiation;
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
% ct_period_initiation = 38 + 1;
% ct_period_refractory = 20 + 1;
% period_initiation = pool_period_initiation( ct_period_initiation );
% period_active = 120 - period_initiation;
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
% ct_period_initiation = 15 + 1;
% ct_period_refractory = 100 + 1;
% period_initiation = pool_period_initiation( ct_period_initiation );
% period_active = 120 - period_initiation;
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
% %
% %
% % ct_period_initiation = 16 + 1;
% % ct_period_refractory = 35 + 1;
% % period_initiation = pool_period_initiation( ct_period_initiation );
% % period_active = 120 - period_initiation;
% % period_refractory = pool_period_refractory( ct_period_refractory );
% %
% % S1 = [];
% % for k = 1 : n_stim
% %     tS = [ zeros( n_stim * pop_size, period_active ) ];
% %     tS( ( k - 1 ) * pop_size + [ 1 : pop_size ], : ) = 1;
% %     S1 = [ S1, tS ];
% % end; clear k tS
% % S2 = circshift( S1, period_initiation, 2);
% % S3 = circshift( S2, period_initiation, 2);
% %
% %
% % figure( 'position', [ 100, 100, 300, 150 ] )
% % subplot( 1, 3, 1 )
% % imagesc( S1 )
% % colormap( gca, [ 1, 1, 1; 0, 0, 0 ] )
% % set( gca, 'ytick', [] )
% % ylabel( 'Assembly' )
% % xlabel( 'Time (ms)' )
% % title( 'Area X' )
% % subplot( 1, 3, 2 )
% % imagesc( S2 )
% % colormap( gca, [ 1, 1, 1; 0, 0, 0 ] )
% % set( gca, 'ytick', [] )
% % xlabel( 'Time (ms)' )
% % title( 'Area Y' )
% % subplot( 1, 3, 3 )
% % imagesc( S3 )
% % colormap( gca, [ 1, 1, 1; 0, 0, 0 ] )
% % set( gca, 'ytick', [] )
% % xlabel( 'Time (ms)' )
% % title( 'Area Z' )
% %

%% Get W1 and W2 for Example test (N = 4)
%
% pool_period_initiation = [ 20 : 60 ];
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
%     period_active = 120 - period_initiation;
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
%     S4 = circshift( S3, period_initiation, 2);
%
%     tic
%     [ W1, W2, W3, W4 ] = learningHebbian4( S1, S2, S3, S4 );
%     toc
%
%     Ws{ ct_period_initiation, 1 } = W1;
%     Ws{ ct_period_initiation, 2 } = W2;
%     Ws{ ct_period_initiation, 3 } = W3;
%     Ws{ ct_period_initiation, 4 } = W4;
%
%     clear period_active S1 S2 S3 S4 W1 W2 W3 W4
%
%     disp( [ num2str( period_initiation ), ' / ', num2str( pool_period_initiation( end ) ), ' parameters.' ] )
%
% end; clear period_initiation ct_period_initiation
%
%
% save( 'W1_W2_W3_W4_example.mat', 'pool_period_initiation', 'Ws', 'pop_size', 'n_stim', 'n_iter_learning' )
%

%% Example test (N = 4)
%
% load( 'W1_W2_W3_W4_example.mat' )
%
% threshold = 0.15;
%
% S = zeros( pop_size * n_stim, n_stim );
% for k = 1 : n_stim
%     S( ( k - 1 ) * pop_size + [ 1 : pop_size ], k ) = 1;
% end; clear k
%
% pool_period_refractory = [ 0 : 120 ];
%
% M = [];
% C = [];
% H = [];
%
% ct_period_initiation = 0;
% for period_initiation = pool_period_initiation
%     ct_period_initiation = ct_period_initiation + 1;
%
%     period_active = 120 - period_initiation;
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
%         S4 = circshift( S3, period_initiation, 2);
%
%         W1 = Ws{ ct_period_initiation, 1 };
%         W2 = Ws{ ct_period_initiation, 2 };
%         W3 = Ws{ ct_period_initiation, 3 };
%         W4 = Ws{ ct_period_initiation, 4 };
%
%         ct_tau = 0;
%         for tau = [ 0 : period_active : size( S1, 2 ) - 1 ]
%             ct_tau = ct_tau + 1;
%
%             tS1 = [ circshift( S1, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS2 = [ circshift( S2, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS3 = [ circshift( S3, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS4 = [ circshift( S4, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%
%             [ Z1, Z2, Z3, Z4 ] = assemblyComputation4( tS1, tS2, tS3, tS4, W1, W2, W3, W4, threshold, period_initiation, period_active, period_refractory );
%
%             M1 = mean( sum( Z1( :, size( S1, 2 ) + 1 :  end ), 1 ), 2 );
%             M2 = mean( sum( Z2( :, size( S2, 2 ) + 1 :  end ), 1 ), 2 );
%             M3 = mean( sum( Z3( :, size( S3, 2 ) + 1 :  end ), 1 ), 2 );
%             M4 = mean( sum( Z4( :, size( S4, 2 ) + 1 :  end ), 1 ), 2 );
%             [ ~, C1, H1 ] = patternMatching( S, Z1( :, size( S1, 2 ) + 1 :  end ) );
%             [ ~, C2, H2 ] = patternMatching( S, Z2( :, size( S2, 2 ) + 1 :  end ) );
%             [ ~, C3, H3 ] = patternMatching( S, Z3( :, size( S3, 2 ) + 1 :  end ) );
%             [ ~, C4, H4 ] = patternMatching( S, Z4( :, size( S4, 2 ) + 1 :  end ) );
%             C1 = C1( H1 > 0 );
%             C2 = C2( H2 > 0 );
%             C3 = C3( H3 > 0 );
%             C4 = C4( H4 > 0 );
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
%             if length( C4 ) >= 0.5 * length( H4 )
%                 t1C = C4( 1 : end - 1 );
%                 t2C = C4( 2 : end );
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
%                 C4 = mean( idxC, 2 ) >= 0.5;
%                 clear t1C t2C idxC
%             else
%                 C4 = false;
%             end
%             H1 = mean( H1, 2 );
%             H2 = mean( H2, 2 );
%             H3 = mean( H3, 2 );
%             H4 = mean( H4, 2 );
%
%             M( ct_period_initiation, ct_period_refractory, ct_tau ) = mean( [ M1, M2, M3, M4 ], 2 );
%             C( ct_period_initiation, ct_period_refractory, ct_tau ) = all( [ C1, C2, C3, C4 ], 2 );
%             H( ct_period_initiation, ct_period_refractory, ct_tau ) = mean( [ H1, H2, H3, H4 ], 2 );
%
%             clear tS1 tS2 tS3 tS4 Z1 Z2 Z3 Z4 M1 M2 M3 M4 C1 C2 C3 C4 H1 H2 H3 H4
%         end; clear tau ct_tau
%
%         clear S1 S2 S3 S4 W1 W2 W3 W4
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
% save( 'Results_example4.mat', 'M', 'C', 'H', 'threshold', 'pool_period_refractory', 'pool_period_initiation', 'Ws', 'pop_size', 'n_stim', 'n_iter_learning' )
%

%% Results for Example test (N = 4)
%
% load( 'Results_example4.mat' )
%
%
% % -------------------------------------------------------------------------
%
% N = 4;
%
% likelihood_f = ( ( 1 - ( transpose( pool_period_initiation ) ./ ( 120 - transpose( pool_period_initiation ) ) ) ) .^ ( N - 1 ) ) .* ( ( N - 1 ) * transpose( pool_period_initiation ) ./ ( 120 - transpose( pool_period_initiation ) ) );
% likelihood_f = repmat( likelihood_f, [ 1, length( pool_period_refractory ) ] );
%
% likelihood_f( repmat( transpose( pool_period_initiation ), [ 1, length( pool_period_refractory ) ] ) <= 0 ) = NaN;
% likelihood_f( bsxfun( @gt, transpose( pool_period_initiation ), ( 120 - pool_period_initiation ) / ( N - 1 ) ) ) = NaN;
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
% P( M > 20 ) = 2;
% P( M >= 1 & M <= 20 & C == 0 ) = 3;
% P( M >= 1 & M <= 20 & C == 1 ) = 4;
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
% ylabel( 'Initiation period (ms)' )
% xlabel( 'Inhibition period (ms)' )
% title( 'States of parameters' )
%

%% Get W1 and W2 for Example test (N = 5)
%
% pool_period_initiation = [ 20 : 60 ];
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
%     period_active = 120 - period_initiation;
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
%     S4 = circshift( S3, period_initiation, 2);
%     S5 = circshift( S4, period_initiation, 2);
%
%     tic
%     [ W1, W2, W3, W4, W5 ] = learningHebbian5( S1, S2, S3, S4, S5 );
%     toc
%
%     Ws{ ct_period_initiation, 1 } = W1;
%     Ws{ ct_period_initiation, 2 } = W2;
%     Ws{ ct_period_initiation, 3 } = W3;
%     Ws{ ct_period_initiation, 4 } = W4;
%     Ws{ ct_period_initiation, 5 } = W5;
%
%     clear period_active S1 S2 S3 S4 S5 W1 W2 W3 W4 W5
%
%     disp( [ num2str( period_initiation ), ' / ', num2str( pool_period_initiation( end ) ), ' parameters.' ] )
%
% end; clear period_initiation ct_period_initiation
%
%
% save( 'W1_W2_W3_W4_W5_example.mat', 'pool_period_initiation', 'Ws', 'pop_size', 'n_stim', 'n_iter_learning' )
%

%% Example test (N = 5)
%
% load( 'W1_W2_W3_W4_W5_example.mat' )
%
% threshold = 0.15;
%
% S = zeros( pop_size * n_stim, n_stim );
% for k = 1 : n_stim
%     S( ( k - 1 ) * pop_size + [ 1 : pop_size ], k ) = 1;
% end; clear k
%
% pool_period_refractory = [ 0 : 120 ];
%
% M = [];
% C = [];
% H = [];
%
% ct_period_initiation = 0;
% for period_initiation = pool_period_initiation
%     ct_period_initiation = ct_period_initiation + 1;
%
%     period_active = 120 - period_initiation;
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
%         S4 = circshift( S3, period_initiation, 2);
%         S5 = circshift( S4, period_initiation, 2);
%
%         W1 = Ws{ ct_period_initiation, 1 };
%         W2 = Ws{ ct_period_initiation, 2 };
%         W3 = Ws{ ct_period_initiation, 3 };
%         W4 = Ws{ ct_period_initiation, 4 };
%         W5 = Ws{ ct_period_initiation, 5 };
%
%         ct_tau = 0;
%         for tau = [ 0 : period_active : size( S1, 2 ) - 1 ]
%             ct_tau = ct_tau + 1;
%
%             tS1 = [ circshift( S1, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS2 = [ circshift( S2, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS3 = [ circshift( S3, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS4 = [ circshift( S4, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS5 = [ circshift( S5, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%
%             [ Z1, Z2, Z3, Z4, Z5 ] = assemblyComputation5( tS1, tS2, tS3, tS4, tS5, W1, W2, W3, W4, W5, threshold, period_initiation, period_active, period_refractory );
%
%             M1 = mean( sum( Z1( :, size( S1, 2 ) + 1 :  end ), 1 ), 2 );
%             M2 = mean( sum( Z2( :, size( S2, 2 ) + 1 :  end ), 1 ), 2 );
%             M3 = mean( sum( Z3( :, size( S3, 2 ) + 1 :  end ), 1 ), 2 );
%             M4 = mean( sum( Z4( :, size( S4, 2 ) + 1 :  end ), 1 ), 2 );
%             M5 = mean( sum( Z5( :, size( S5, 2 ) + 1 :  end ), 1 ), 2 );
%             [ ~, C1, H1 ] = patternMatching( S, Z1( :, size( S1, 2 ) + 1 :  end ) );
%             [ ~, C2, H2 ] = patternMatching( S, Z2( :, size( S2, 2 ) + 1 :  end ) );
%             [ ~, C3, H3 ] = patternMatching( S, Z3( :, size( S3, 2 ) + 1 :  end ) );
%             [ ~, C4, H4 ] = patternMatching( S, Z4( :, size( S4, 2 ) + 1 :  end ) );
%             [ ~, C5, H5 ] = patternMatching( S, Z5( :, size( S5, 2 ) + 1 :  end ) );
%             C1 = C1( H1 > 0 );
%             C2 = C2( H2 > 0 );
%             C3 = C3( H3 > 0 );
%             C4 = C4( H4 > 0 );
%             C5 = C5( H5 > 0 );
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
%             if length( C4 ) >= 0.5 * length( H4 )
%                 t1C = C4( 1 : end - 1 );
%                 t2C = C4( 2 : end );
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
%                 C4 = mean( idxC, 2 ) >= 0.5;
%                 clear t1C t2C idxC
%             else
%                 C4 = false;
%             end
%             if length( C5 ) >= 0.5 * length( H4 )
%                 t1C = C5( 1 : end - 1 );
%                 t2C = C5( 2 : end );
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
%                 C5 = mean( idxC, 2 ) >= 0.5;
%                 clear t1C t2C idxC
%             else
%                 C5 = false;
%             end
%             H1 = mean( H1, 2 );
%             H2 = mean( H2, 2 );
%             H3 = mean( H3, 2 );
%             H4 = mean( H4, 2 );
%             H5 = mean( H5, 2 );
%
%             M( ct_period_initiation, ct_period_refractory, ct_tau ) = mean( [ M1, M2, M3, M4, M5 ], 2 );
%             C( ct_period_initiation, ct_period_refractory, ct_tau ) = all( [ C1, C2, C3, C4, C5 ], 2 );
%             H( ct_period_initiation, ct_period_refractory, ct_tau ) = mean( [ H1, H2, H3, H4, H5 ], 2 );
%
%             clear tS1 tS2 tS3 tS4 tS5 Z1 Z2 Z3 Z4 Z5 M1 M2 M3 M4 M5 C1 C2 C3 C4 C5 H1 H2 H3 H4 H5
%         end; clear tau ct_tau
%
%         clear S1 S2 S3 S4 S5 W1 W2 W3 W4 W5
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
% save( 'Results_example5.mat', 'M', 'C', 'H', 'threshold', 'pool_period_refractory', 'pool_period_initiation', 'Ws', 'pop_size', 'n_stim', 'n_iter_learning' )
%

%% Results for Example test (N = 5)
%
% load( 'Results_example5.mat' )
%
%
% % -------------------------------------------------------------------------
%
% N = 5;
%
% likelihood_f = ( ( 1 - ( transpose( pool_period_initiation ) ./ ( 120 - transpose( pool_period_initiation ) ) ) ) .^ ( N - 1 ) ) .* ( ( N - 1 ) * transpose( pool_period_initiation ) ./ ( 120 - transpose( pool_period_initiation ) ) );
% likelihood_f = repmat( likelihood_f, [ 1, length( pool_period_refractory ) ] );
%
% likelihood_f( repmat( transpose( pool_period_initiation ), [ 1, length( pool_period_refractory ) ] ) <= 0 ) = NaN;
% likelihood_f( bsxfun( @gt, transpose( pool_period_initiation ), ( 120 - pool_period_initiation ) / ( N - 1 ) ) ) = NaN;
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
% P( M > 20 ) = 2;
% P( M >= 1 & M <= 20 & C == 0 ) = 3;
% P( M >= 1 & M <= 20 & C == 1 ) = 4;
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
% ylabel( 'Initiation period (ms)' )
% xlabel( 'Inhibition period (ms)' )
% title( 'States of parameters' )
%

%% Get W1 and W2 for Example test (N = 6)
%
% pool_period_initiation = [ 20 : 60 ];
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
%     period_active = 120 - period_initiation;
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
%     S4 = circshift( S3, period_initiation, 2);
%     S5 = circshift( S4, period_initiation, 2);
%     S6 = circshift( S5, period_initiation, 2);
%
%     tic
%     [ W1, W2, W3, W4, W5, W6 ] = learningHebbian6( S1, S2, S3, S4, S5, S6 );
%     toc
%
%     Ws{ ct_period_initiation, 1 } = W1;
%     Ws{ ct_period_initiation, 2 } = W2;
%     Ws{ ct_period_initiation, 3 } = W3;
%     Ws{ ct_period_initiation, 4 } = W4;
%     Ws{ ct_period_initiation, 5 } = W5;
%     Ws{ ct_period_initiation, 6 } = W6;
%
%     clear period_active S1 S2 S3 S4 S5 S6 W1 W2 W3 W4 W5 W6
%
%     disp( [ num2str( period_initiation ), ' / ', num2str( pool_period_initiation( end ) ), ' parameters.' ] )
%
% end; clear period_initiation ct_period_initiation
%
%
% save( 'W1_W2_W3_W4_W5_W6_example.mat', 'pool_period_initiation', 'Ws', 'pop_size', 'n_stim', 'n_iter_learning' )
%

%% Example test (N = 6)
%
% load( 'W1_W2_W3_W4_W5_W6_example.mat' )
%
% threshold = 0.15;
%
% S = zeros( pop_size * n_stim, n_stim );
% for k = 1 : n_stim
%     S( ( k - 1 ) * pop_size + [ 1 : pop_size ], k ) = 1;
% end; clear k
%
% pool_period_refractory = [ 0 : 120 ];
%
% M = [];
% C = [];
% H = [];
%
% ct_period_initiation = 0;
% for period_initiation = pool_period_initiation
%     ct_period_initiation = ct_period_initiation + 1;
%
%     period_active = 120 - period_initiation;
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
%         S4 = circshift( S3, period_initiation, 2);
%         S5 = circshift( S4, period_initiation, 2);
%         S6 = circshift( S5, period_initiation, 2);
%
%         W1 = Ws{ ct_period_initiation, 1 };
%         W2 = Ws{ ct_period_initiation, 2 };
%         W3 = Ws{ ct_period_initiation, 3 };
%         W4 = Ws{ ct_period_initiation, 4 };
%         W5 = Ws{ ct_period_initiation, 5 };
%         W6 = Ws{ ct_period_initiation, 6 };
%
%         ct_tau = 0;
%         for tau = [ 0 : period_active : size( S1, 2 ) - 1 ]
%             ct_tau = ct_tau + 1;
%
%             tS1 = [ circshift( S1, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS2 = [ circshift( S2, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS3 = [ circshift( S3, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS4 = [ circshift( S4, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS5 = [ circshift( S5, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS6 = [ circshift( S6, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%
%             [ Z1, Z2, Z3, Z4, Z5, Z6 ] = assemblyComputation6( tS1, tS2, tS3, tS4, tS5, tS6, W1, W2, W3, W4, W5, W6, threshold, period_initiation, period_active, period_refractory );
%
%             M1 = mean( sum( Z1( :, size( S1, 2 ) + 1 :  end ), 1 ), 2 );
%             M2 = mean( sum( Z2( :, size( S2, 2 ) + 1 :  end ), 1 ), 2 );
%             M3 = mean( sum( Z3( :, size( S3, 2 ) + 1 :  end ), 1 ), 2 );
%             M4 = mean( sum( Z4( :, size( S4, 2 ) + 1 :  end ), 1 ), 2 );
%             M5 = mean( sum( Z5( :, size( S5, 2 ) + 1 :  end ), 1 ), 2 );
%             M6 = mean( sum( Z6( :, size( S6, 2 ) + 1 :  end ), 1 ), 2 );
%             [ ~, C1, H1 ] = patternMatching( S, Z1( :, size( S1, 2 ) + 1 :  end ) );
%             [ ~, C2, H2 ] = patternMatching( S, Z2( :, size( S2, 2 ) + 1 :  end ) );
%             [ ~, C3, H3 ] = patternMatching( S, Z3( :, size( S3, 2 ) + 1 :  end ) );
%             [ ~, C4, H4 ] = patternMatching( S, Z4( :, size( S4, 2 ) + 1 :  end ) );
%             [ ~, C5, H5 ] = patternMatching( S, Z5( :, size( S5, 2 ) + 1 :  end ) );
%             [ ~, C6, H6 ] = patternMatching( S, Z6( :, size( S6, 2 ) + 1 :  end ) );
%             C1 = C1( H1 > 0 );
%             C2 = C2( H2 > 0 );
%             C3 = C3( H3 > 0 );
%             C4 = C4( H4 > 0 );
%             C5 = C5( H5 > 0 );
%             C6 = C6( H6 > 0 );
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
%             if length( C4 ) >= 0.5 * length( H4 )
%                 t1C = C4( 1 : end - 1 );
%                 t2C = C4( 2 : end );
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
%                 C4 = mean( idxC, 2 ) >= 0.5;
%                 clear t1C t2C idxC
%             else
%                 C4 = false;
%             end
%             if length( C5 ) >= 0.5 * length( H4 )
%                 t1C = C5( 1 : end - 1 );
%                 t2C = C5( 2 : end );
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
%                 C5 = mean( idxC, 2 ) >= 0.5;
%                 clear t1C t2C idxC
%             else
%                 C5 = false;
%             end
%             if length( C6 ) >= 0.5 * length( H4 )
%                 t1C = C6( 1 : end - 1 );
%                 t2C = C6( 2 : end );
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
%                 C6 = mean( idxC, 2 ) >= 0.5;
%                 clear t1C t2C idxC
%             else
%                 C6 = false;
%             end
%             H1 = mean( H1, 2 );
%             H2 = mean( H2, 2 );
%             H3 = mean( H3, 2 );
%             H4 = mean( H4, 2 );
%             H5 = mean( H5, 2 );
%             H6 = mean( H6, 2 );
%
%             M( ct_period_initiation, ct_period_refractory, ct_tau ) = mean( [ M1, M2, M3, M4, M5, M6 ], 2 );
%             C( ct_period_initiation, ct_period_refractory, ct_tau ) = all( [ C1, C2, C3, C4, C5, C6 ], 2 );
%             H( ct_period_initiation, ct_period_refractory, ct_tau ) = mean( [ H1, H2, H3, H4, H5, H6 ], 2 );
%
%             clear tS1 tS2 tS3 tS4 tS5 tS6 Z1 Z2 Z3 Z4 Z5 Z6 M1 M2 M3 M4 M5 M6 C1 C2 C3 C4 C5 C6 H1 H2 H3 H4 H5 H6
%         end; clear tau ct_tau
%
%         clear S1 S2 S3 S4 S5 S6 W1 W2 W3 W4 W5 W6
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
% save( 'Results_example6.mat', 'M', 'C', 'H', 'threshold', 'pool_period_refractory', 'pool_period_initiation', 'Ws', 'pop_size', 'n_stim', 'n_iter_learning' )
%

%% Results for Example test (N = 6)
%
% load( 'Results_example6.mat' )
%
%
% % -------------------------------------------------------------------------
%
% N = 6;
%
% likelihood_f = ( ( 1 - ( transpose( pool_period_initiation ) ./ ( 120 - transpose( pool_period_initiation ) ) ) ) .^ ( N - 1 ) ) .* ( ( N - 1 ) * transpose( pool_period_initiation ) ./ ( 120 - transpose( pool_period_initiation ) ) );
% likelihood_f = repmat( likelihood_f, [ 1, length( pool_period_refractory ) ] );
%
% likelihood_f( repmat( transpose( pool_period_initiation ), [ 1, length( pool_period_refractory ) ] ) <= 0 ) = NaN;
% likelihood_f( bsxfun( @gt, transpose( pool_period_initiation ), ( 120 - pool_period_initiation ) / ( N - 1 ) ) ) = NaN;
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
% P( M > 20 ) = 2;
% P( M >= 1 & M <= 20 & C == 0 ) = 3;
% P( M >= 1 & M <= 20 & C == 1 ) = 4;
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
% ylabel( 'Initiation period (ms)' )
% xlabel( 'Inhibition period (ms)' )
% title( 'States of parameters' )
%

%% Get W1 and W2 for Example test (N = 7)
%
% pool_period_initiation = [ 20 : 60 ];
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
%     period_active = 120 - period_initiation;
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
%     S4 = circshift( S3, period_initiation, 2);
%     S5 = circshift( S4, period_initiation, 2);
%     S6 = circshift( S5, period_initiation, 2);
%     S7 = circshift( S6, period_initiation, 2);
%
%     tic
%     [ W1, W2, W3, W4, W5, W6, W7 ] = learningHebbian7( S1, S2, S3, S4, S5, S6, S7 );
%     toc
%
%     Ws{ ct_period_initiation, 1 } = W1;
%     Ws{ ct_period_initiation, 2 } = W2;
%     Ws{ ct_period_initiation, 3 } = W3;
%     Ws{ ct_period_initiation, 4 } = W4;
%     Ws{ ct_period_initiation, 5 } = W5;
%     Ws{ ct_period_initiation, 6 } = W6;
%     Ws{ ct_period_initiation, 7 } = W7;
%
%     clear period_active S1 S2 S3 S4 S5 S6 S7 W1 W2 W3 W4 W5 W6 W7
%
%     disp( [ num2str( period_initiation ), ' / ', num2str( pool_period_initiation( end ) ), ' parameters.' ] )
%
% end; clear period_initiation ct_period_initiation
%
%
% save( 'W1_W2_W3_W4_W5_W6_W7_example.mat', 'pool_period_initiation', 'Ws', 'pop_size', 'n_stim', 'n_iter_learning' )
%

%% Example test (N = 7)
%
% load( 'W1_W2_W3_W4_W5_W6_W7_example.mat' )
%
% threshold = 0.15;
%
% S = zeros( pop_size * n_stim, n_stim );
% for k = 1 : n_stim
%     S( ( k - 1 ) * pop_size + [ 1 : pop_size ], k ) = 1;
% end; clear k
%
% pool_period_refractory = [ 0 : 120 ];
%
% M = [];
% C = [];
% H = [];
%
% ct_period_initiation = 0;
% for period_initiation = pool_period_initiation
%     ct_period_initiation = ct_period_initiation + 1;
%
%     period_active = 120 - period_initiation;
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
%         S4 = circshift( S3, period_initiation, 2);
%         S5 = circshift( S4, period_initiation, 2);
%         S6 = circshift( S5, period_initiation, 2);
%         S7 = circshift( S6, period_initiation, 2);
%
%         W1 = Ws{ ct_period_initiation, 1 };
%         W2 = Ws{ ct_period_initiation, 2 };
%         W3 = Ws{ ct_period_initiation, 3 };
%         W4 = Ws{ ct_period_initiation, 4 };
%         W5 = Ws{ ct_period_initiation, 5 };
%         W6 = Ws{ ct_period_initiation, 6 };
%         W7 = Ws{ ct_period_initiation, 7 };
%
%         ct_tau = 0;
%         for tau = [ 0 : period_active : size( S1, 2 ) - 1 ]
%             ct_tau = ct_tau + 1;
%
%             tS1 = [ circshift( S1, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS2 = [ circshift( S2, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS3 = [ circshift( S3, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS4 = [ circshift( S4, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS5 = [ circshift( S5, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS6 = [ circshift( S6, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%             tS7 = [ circshift( S7, tau, 2 ), NaN( size( S1, 1 ), 3000 ) ];
%
%             [ Z1, Z2, Z3, Z4, Z5, Z6, Z7 ] = assemblyComputation7( tS1, tS2, tS3, tS4, tS5, tS6, tS7, W1, W2, W3, W4, W5, W6, W7, threshold, period_initiation, period_active, period_refractory );
%
%             M1 = mean( sum( Z1( :, size( S1, 2 ) + 1 :  end ), 1 ), 2 );
%             M2 = mean( sum( Z2( :, size( S2, 2 ) + 1 :  end ), 1 ), 2 );
%             M3 = mean( sum( Z3( :, size( S3, 2 ) + 1 :  end ), 1 ), 2 );
%             M4 = mean( sum( Z4( :, size( S4, 2 ) + 1 :  end ), 1 ), 2 );
%             M5 = mean( sum( Z5( :, size( S5, 2 ) + 1 :  end ), 1 ), 2 );
%             M6 = mean( sum( Z6( :, size( S6, 2 ) + 1 :  end ), 1 ), 2 );
%             M7 = mean( sum( Z7( :, size( S7, 2 ) + 1 :  end ), 1 ), 2 );
%             [ ~, C1, H1 ] = patternMatching( S, Z1( :, size( S1, 2 ) + 1 :  end ) );
%             [ ~, C2, H2 ] = patternMatching( S, Z2( :, size( S2, 2 ) + 1 :  end ) );
%             [ ~, C3, H3 ] = patternMatching( S, Z3( :, size( S3, 2 ) + 1 :  end ) );
%             [ ~, C4, H4 ] = patternMatching( S, Z4( :, size( S4, 2 ) + 1 :  end ) );
%             [ ~, C5, H5 ] = patternMatching( S, Z5( :, size( S5, 2 ) + 1 :  end ) );
%             [ ~, C6, H6 ] = patternMatching( S, Z6( :, size( S6, 2 ) + 1 :  end ) );
%             [ ~, C7, H7 ] = patternMatching( S, Z7( :, size( S7, 2 ) + 1 :  end ) );
%             C1 = C1( H1 > 0 );
%             C2 = C2( H2 > 0 );
%             C3 = C3( H3 > 0 );
%             C4 = C4( H4 > 0 );
%             C5 = C5( H5 > 0 );
%             C6 = C6( H6 > 0 );
%             C7 = C7( H7 > 0 );
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
%             if length( C4 ) >= 0.5 * length( H4 )
%                 t1C = C4( 1 : end - 1 );
%                 t2C = C4( 2 : end );
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
%                 C4 = mean( idxC, 2 ) >= 0.5;
%                 clear t1C t2C idxC
%             else
%                 C4 = false;
%             end
%             if length( C5 ) >= 0.5 * length( H4 )
%                 t1C = C5( 1 : end - 1 );
%                 t2C = C5( 2 : end );
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
%                 C5 = mean( idxC, 2 ) >= 0.5;
%                 clear t1C t2C idxC
%             else
%                 C5 = false;
%             end
%             if length( C6 ) >= 0.5 * length( H4 )
%                 t1C = C6( 1 : end - 1 );
%                 t2C = C6( 2 : end );
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
%                 C6 = mean( idxC, 2 ) >= 0.5;
%                 clear t1C t2C idxC
%             else
%                 C6 = false;
%             end
%             if length( C7 ) >= 0.5 * length( H4 )
%                 t1C = C7( 1 : end - 1 );
%                 t2C = C7( 2 : end );
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
%                 C7 = mean( idxC, 2 ) >= 0.5;
%                 clear t1C t2C idxC
%             else
%                 C7 = false;
%             end
%             H1 = mean( H1, 2 );
%             H2 = mean( H2, 2 );
%             H3 = mean( H3, 2 );
%             H4 = mean( H4, 2 );
%             H5 = mean( H5, 2 );
%             H6 = mean( H6, 2 );
%             H7 = mean( H7, 2 );
%
%             M( ct_period_initiation, ct_period_refractory, ct_tau ) = mean( [ M1, M2, M3, M4, M5, M6, M7 ], 2 );
%             C( ct_period_initiation, ct_period_refractory, ct_tau ) = all( [ C1, C2, C3, C4, C5, C6, C7 ], 2 );
%             H( ct_period_initiation, ct_period_refractory, ct_tau ) = mean( [ H1, H2, H3, H4, H5, H6, H7 ], 2 );
%
%             clear tS1 tS2 tS3 tS4 tS5 tS6 tS7 Z1 Z2 Z3 Z4 Z5 Z6 Z7 M1 M2 M3 M4 M5 M6 M7 C1 C2 C3 C4 C5 C6 C7 H1 H2 H3 H4 H5 H6 H7
%         end; clear tau ct_tau
%
%         clear S1 S2 S3 S4 S5 S6 S7 W1 W2 W3 W4 W5 W6 W7
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
% save( 'Results_example7.mat', 'M', 'C', 'H', 'threshold', 'pool_period_refractory', 'pool_period_initiation', 'Ws', 'pop_size', 'n_stim', 'n_iter_learning' )
%

%% SNN_Simul2_0
%
% clear vars
%
% vars.dt = 0.5;% time resolution (ms)
% vars.t_memory = 250;% time memory (ms)
%
% vars.tau_E = 20;% E neuron resting membrane time constant (ms)
% vars.tau_I = 20;% I neuron resting membrane time constant (ms)
% vars.E_L_E = -70;% E neuron resting potential (mV)
% vars.E_L_I = -62;% I neuron resting potential (mV)
% vars.Delta_T_E = 2;% E neuron EIF slope factor (mV)
% vars.Delta_T_I = 2;% I neuron EIF slope factor (mV)
% vars.C = 300;% Capacitance (pF)
% vars.E_E = 0;% E reverse potential (mV)
% vars.E_I = -75;% I reverse potential (mV)
% vars.V_T = -52;% threshold potential (mV)
% vars.A_T = 10;% Post spike threshold potential increase (mV)
% vars.tau_T = 30;% Adaptive threshold time scale (ms)
% vars.V_re = -60;% Reset potential (mV)
% vars.tau_abs = 1;% Absolute refractory period (ms)
% vars.a_w = 4;% Subthreshold adaptation (nS)
% vars.b_w = 0.805;% Spike-triggered adaptation (pA)
% vars.tau_w = 150;% Spike-triggered adaptation time scale (ms)
% vars.V_ap = 20;% action potential threshold potention (mV)
% vars.V_lb = -75;% lower bound potential (mV) := vars.E_I
%
% vars.N_E = [ 500, 500 ];% the number of excitatory neurons
% vars.N_I = [ 125, 125 ];% the number of excitatory neurons
% % vars.N_E = [ 1000, 1000 ];% the number of excitatory neurons
% % vars.N_I = [ 250, 250 ];% the number of excitatory neurons
% % vars.p_C = 2.5 * 0.2;% connection probability via chemical synapses
% vars.p_C = 0.2;% connection probability via chemical synapses
% vars.scaleWeight = 0.14;% parameter for re-scaling of synaptic weights
% % vars.scaleWeight = 0.12;% parameter for re-scaling of synaptic weights
% % vars.scaleWeight = 0.15;% parameter for re-scaling of synaptic weights
% % vars.scaleWeight = 0.2;% parameter for re-scaling of synaptic weights
% % vars.p_E = 0.2;% connection probability via electrical synapses
% % vars.tau_r_E = 1;% rise time for E synapses (ms)
% % vars.tau_d_E = 6;% decay time for E synapses (ms)
% % vars.tau_r_I = 0.5;% rise time for I synapses (ms)
% % vars.tau_d_I = 2;% decay time for I synapses (ms)
% vars.r_ext_EE = 12 * 500;% Rate of external input to E neurons (Hz)
% vars.r_ext_IE = 0 * 500;% Rate of external input to I neurons (Hz)
% vars.J_EE_min = 1.78;% minimum E to E synaptic weight (pF)
% vars.J_EE_max = 21.4;% maximum E to E synaptic weight (pF)
% vars.J_EE_O = 2.76;% initial E to E synaptic weight (pF)
% vars.J_EE_O1 = 2.76 * 8;% initial E to E synaptic weight (pF)
% vars.J_EI_min = 48.7;% minimum I to E synaptic weight (pF)
% vars.J_EI_max = 243;% maximum I to E synaptic weight (pF)
% vars.J_EI_O = 3 * 48.7;% initial I to E synaptic weight (pF)
% vars.J_IE_val = 1.27;% value of E to I synaptic weight (pF)
% vars.J_II_val = 16.2;% value of I to I synaptic weight (pF)
% % vars.J_GII_val = 0;% value of I to I gap junction weight (pF)
% vars.J_ext_EE = vars.J_EE_min;
% vars.J_ext_IE = vars.J_IE_val;
% vars.period_scaling = 20;% period of synaptic scaling (ms)
%
% vars.A_LTD = 0.0008;% Long-term depression (LTD) strength (pA mV^-1)
% vars.A_LTP = 0.0014;% Long-term potentiation (LTP) strength (pA mV^-2)
% vars.theta_LTD = -70;% Threshold to recruit LTD (mV)
% vars.theta_LTP = -49;% Threshold to recruit LTP (mV)
% vars.tau_u = 10;% Time constant of low-pass filtered membrane voltage (for LTD) (ms)
% vars.tau_v = 7;% Time constant of low-pass filtered membrane voltage (for LTP) (ms)
% vars.tau_x = 15;% Time constant low-pass filtered spike train (for LTP) (ms)
%
% vars.tau_y = 20;% Time constant of low-pass filtered spike train (ms)
% vars.eta = 1;% Synaptic plasticity learning rate (pA)
% vars.r_O = 3;% Target firing rate (Hz)
%
% % vars.triplet_tau_plus = 16.8;% parameter of triplet plasticity (ms)
% % vars.triplet_tau_minus = 33.7;% parameter of triplet plasticity (ms)
% % vars.triplet_tau_x = 101;% parameter of triplet plasticity (ms)
% % vars.triplet_tau_y = 125;% parameter of triplet plasticity (ms)
% % vars.triplet_A_2_plus = 7.5 * 10 ^ (-10);% parameter of triplet plasticity
% % vars.triplet_A_3_plus = 9.3 * 10 ^ (-3);% parameter of triplet plasticity
% % vars.triplet_A_2_minus = 7 * 10 ^ (-3);% parameter of triplet plasticity
% % vars.triplet_A_3_minus = 2.3 * 10 ^ (-4);% parameter of triplet plasticity
%
% vars.V_E = [];% voltages of E (mV)
% vars.V_I = [];% voltages of I (mV)
% vars.Vu_E = [];% low-pass filtered valtage of E (mV)
% vars.Vv_E = [];% low-pass filtered valtage of E (mV)
%
% vars.s_E = [];% spike train of E
% vars.sx_E = [];% low-pass filtered spike train of E
% vars.sy_E = [];% low-pass filtered spike train of E
% vars.s_I = [];% spike train of I
% vars.sy_I = [];% low-pass filtered spike train of I
%
% vars.V_T_E = [];% thereshold of E
% vars.ww_E = [];% adaptation current of E
%
% vars.s_ext_EE = [];% external stimulation to E
% vars.s_ext_IE = [];% external stimulation to I
%
% vars.ct_abs_E = [];% count of absolute refractory period for E neurons
% vars.ct_abs_I = [];% count of absolute refractory period for I neurons
%
% % vars.triplet_r_1 = [];% variable of triplet plasticity
% % vars.triplet_r_2 = [];% variable of triplet plasticity
% % vars.triplet_o_1 = [];% variable of triplet plasticity
% % vars.triplet_o_2 = [];% variable of triplet plasticity
%
% vars.connScale = [ 1, 1, 7 ];% [ orthodromic, antidromic, within area ]
% vars.connects = zeros( length( vars.N_E ), length( vars.N_E ) );% vars.connects( j, i ): from area i to area j.
% % vars.connects( 2, 1 ) = vars.connScale( 1 ) * vars.N_E( 1 );
% % vars.connects( 1, 2 ) = vars.connScale( 2 ) * vars.N_E( 2 );
% % vars.connects( 1, 1 ) = vars.connScale( 3 ) * vars.N_E( 1 );
% % vars.connects( 2, 2 ) = vars.connScale( 3 ) * vars.N_E( 2 );
% % vars.connects = vars.connects ./ sum( vars.connects, 2 );
% % vars.connects( isnan( vars.connects ) ) = 0;
% % vars.connects( vars.connects > 0 ) = 1;
% vars.connects( 2, 1 ) = vars.connScale( 1 );
% vars.connects( 1, 2 ) = vars.connScale( 2 );
% vars.connects( 1, 1 ) = vars.connScale( 3 );
% vars.connects( 2, 2 ) = vars.connScale( 3 );
%
% vars.tau_r_E = 1 * 2;% rise time for E synapses (ms)
% vars.tau_d_E = 6 * 2;% decay time for E synapses (ms)
% vars.tau_r_I = 0.5 * 10;% rise time for I synapses (ms)
% vars.tau_d_I = 2 * 10;% decay time for I synapses (ms)
% vars.nChannels = 10;
% vars = SNN_initiation0( vars );
%
%
% save( 'SNN_Simul2_0.mat', 'vars' )
% clear vars
%
%
% % dt = vars.dt;
% % t_memory = vars.t_memory;
% % t_memory_L = t_memory / dt;
% % tau_r_E = vars.tau_r_E;
% % tau_d_E = vars.tau_d_E;
% % tau_r_I = vars.tau_r_I;
% % tau_d_I = vars.tau_d_I;
% % t_domain_memory = 0 : dt : t_memory - dt;
% % F_E = ( 1 / ( tau_d_E - tau_r_E ) ) * ( exp( -t_domain_memory / tau_d_E ) - exp( -t_domain_memory / tau_r_E ) );
% % F_I = ( 1 / ( tau_d_I - tau_r_I ) ) * ( exp( -t_domain_memory / tau_d_I ) - exp( -t_domain_memory / tau_r_I ) );
% % figure; hold on; plot(F_E); plot(F_I)
%

%% SNN_Simul2_1
%
% transMatU = zeros( 20, 1000 );
% ct_u = 0;
% for u = 1 : 20
%     transMatU( u, ct_u + [ 1 : 50 ] ) = 1;
%     ct_u = ct_u + 50;
% end; clear u ct_u
%
% transMatT = zeros( 9600, 480 );
% ct_t = 0;
% for t = 1 : 480
%     transMatT( ct_t + [ 1 : 20 ], t ) = 1;
%     ct_t = ct_t + 20;
% end; clear t ct_t
%
% n_stim = 10;
% S = eye( n_stim );
%
%
% pool_E = [ 7 : 2 : 25 ];
% pool_I = [ 3 : 3 : 30 ];
%
% M = [];
% C = [];
% H = [];
% s_E_all = {};
% s_E_all_trans = {};
%
% ct_E = 0;
% for mE = pool_E
%     ct_E = ct_E + 1;
%     ct_I = 0;
%     for mI = pool_I
%         ct_I = ct_I + 1;
%
%         load( 'SNN_Simul2_0.mat' )
%
%         vars.tau_r_E = 1 * mE;% rise time for E synapses (ms)
%         vars.tau_d_E = 6 * mE * 0.3;% decay time for E synapses (ms)
%         vars.tau_r_I = 0.5 * mI;% rise time for I synapses (ms)
%         vars.tau_d_I = 2 * mI;% decay time for I synapses (ms)
%
%         vars.mode = 0.5;% the current balance between synaptic weights ( 1:= lateral, 0:= recurrent ), \in [ 0, 1 ].
%         vars.infoStims = [ 4800, 120, 40, 10, 1 ];% recording duration (ms), stimulation duration (ms), stimulation delay (ms), the number of assemblies, the number of stimulations ];
%         vars.inputs = [ 1, 2 ];% in area
%         vars = SNN_stimulation( vars, 0 );
%         [ vars, records ] = SNN_simulation( vars, 'plasticity_off', 'records_on' );
%
%         s_E = records.s_E;
%         s_E_all{ ct_E, ct_I } = s_E;
%         s_E = ( transMatU * s_E ) * transMatT;
%         s_E_all_trans{ ct_E, ct_I } = s_E;
%
%         % figure
%         % imagesc( vars.ext_stim_E ); colorbar
%         % figure
%         % imagesc( records.s_E ); colorbar
%         % figure
%         % imagesc( s_E ); colorbar
%
%
%         Z1 = s_E( 1 : 10, : );
%         Z2 = s_E( 11 : 20, : );
%
%         M1 = mean( sum( Z1, 1 ), 2 );
%         M2 = mean( sum( Z2, 1 ), 2 );
%         [ ~, C1, H1 ] = patternMatching( S, Z1 );
%         [ ~, C2, H2 ] = patternMatching( S, Z2 );
%         C1 = C1( H1 > 0 );
%         C2 = C2( H2 > 0 );
%         if length( C1 ) >= 0.5 * length( H1 )
%             t1C = C1( 1 : end - 1 );
%             t2C = C1( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C1 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C1 = false;
%         end
%         if length( C2 ) >= 0.5 * length( H2 )
%             t1C = C2( 1 : end - 1 );
%             t2C = C2( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C2 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C2 = false;
%         end
%         H1 = mean( H1, 2 );
%         H2 = mean( H2, 2 );
%
%         M( ct_E, ct_I ) = mean( [ M1, M2 ], 2 );
%         C( ct_E, ct_I ) = all( [ C1, C2 ], 2 );
%         H( ct_E, ct_I ) = mean( [ H1, H2 ], 2 );
%
%         disp( [ 'E: ', num2str( ct_E ), ' / I: ', num2str( ct_I ) ] )
%     end; clear mI ct_I
% end; clear mE ct_E
%
% save( 'SNN_Simul2_1.mat', 'M', 'C', 'H', 's_E_all', 's_E_all_trans', 'pool_E', 'pool_I', '-v7.3' )
%

%% Result, SNN_Simul2
%
% load( 'SNN_Simul2_0.mat' )
% load( 'SNN_Simul2_1.mat' )
%
%
% figure( 'position', [ 100, 100, 300, 220 ] )
% imagesc( vars.J_EE, [ 0, 20 ] )
% colorbar% (pF)
% set( gca, 'xtick', [], 'ytick', [] )
%
% % ct_f = 0;
% % figure
% % for ct_E = 1 : 10
% %     for ct_I = 1 : 10
% %         ct_f = ct_f + 1;
% %         subplot( 10, 10, ct_f )
% %         imagesc( s_E_all_trans{ ct_E, ct_I } )
% %     end; clear ct_I
% % end; clear ct_E
%
%
% transMatU = zeros( 20, 1000 );
% ct_u = 0;
% for u = 1 : 20
%     transMatU( u, ct_u + [ 1 : 50 ] ) = 1;
%     ct_u = ct_u + 50;
% end; clear u ct_u
%
% spkT1 = [];% first spike timing
% for ct_E = 1 : length( pool_E )
%     for ct_I = 1 : length( pool_I )
%         s_E = s_E_all{ ct_E, ct_I };
%         s_E = ( transMatU * s_E );
%         spkT1( ct_E, ct_I ) = find( s_E( 1, : ) > 25, 1 );
%         clear s_E
%     end; clear ct_I
% end; clear ct_E
% spkT1 = mean( spkT1, 2 ) / 2;% (ms)
%
%
% colors = [ 1, 1, 1; 0.8, 0.8, 0.3 ];
%
% figure( 'position', [ 100, 100, 300, 200 ] )
% imagesc( 2 * ( 0.5 + 2 ) * pool_I, 1 * pool_E, C, [ -0.5, 1.5 ] )
% colormap( gca, colors )
% colorbar( 'ticks', [ 0, 1 ], 'ticklabels', { '', '' } )
% axis xy
% hold on
% ylabel( 'Excitatory rise param. (ms)' )
% xlabel( 'Approx. inhibition dur. (ms)' )
% title( 'States of parameters' )
%
% lw = 2;
% figure( 'position', [ 100, 100, 200, 200 ] )
% plot( pool_E, spkT1, '-', 'linewidth', lw )
% xlabel( 'Excitatory rise param. (ms)' )
% ylabel( 'Actual initiation dur. (ms)' )
%
% figure( 'position', [ 100, 100, 300, 200 ] )
% imagesc( 2 * ( 0.5 + 2 ) * pool_I, spkT1, C, [ -0.5, 1.5 ] )
% colormap( gca, colors )
% colorbar( 'ticks', [ 0, 1 ], 'ticklabels', { '', '' } )
% axis xy
% hold on
% ylabel( 'Actual initiation dur. (ms)' )
% xlabel( 'Approx. inhibition dur. (ms)' )
% title( 'States of parameters' )
%

%% SNN_Simul3_0
%
% clear vars
%
% vars.dt = 0.5;% time resolution (ms)
% vars.t_memory = 250;% time memory (ms)
%
% vars.tau_E = 20;% E neuron resting membrane time constant (ms)
% vars.tau_I = 20;% I neuron resting membrane time constant (ms)
% vars.E_L_E = -70;% E neuron resting potential (mV)
% vars.E_L_I = -62;% I neuron resting potential (mV)
% vars.Delta_T_E = 2;% E neuron EIF slope factor (mV)
% vars.Delta_T_I = 2;% I neuron EIF slope factor (mV)
% vars.C = 300;% Capacitance (pF)
% vars.E_E = 0;% E reverse potential (mV)
% vars.E_I = -75;% I reverse potential (mV)
% vars.V_T = -52;% threshold potential (mV)
% vars.A_T = 10;% Post spike threshold potential increase (mV)
% vars.tau_T = 30;% Adaptive threshold time scale (ms)
% vars.V_re = -60;% Reset potential (mV)
% vars.tau_abs = 1;% Absolute refractory period (ms)
% vars.a_w = 4;% Subthreshold adaptation (nS)
% vars.b_w = 0.805;% Spike-triggered adaptation (pA)
% vars.tau_w = 150;% Spike-triggered adaptation time scale (ms)
% vars.V_ap = 20;% action potential threshold potention (mV)
% vars.V_lb = -75;% lower bound potential (mV) := vars.E_I
%
% vars.N_E = [ 500, 500, 500 ];% the number of excitatory neurons
% vars.N_I = [ 125, 125, 125 ];% the number of excitatory neurons
% % vars.N_E = [ 1000, 1000 ];% the number of excitatory neurons
% % vars.N_I = [ 250, 250 ];% the number of excitatory neurons
% % vars.p_C = 2.5 * 0.2;% connection probability via chemical synapses
% vars.p_C = 0.2;% connection probability via chemical synapses
% vars.scaleWeight = 0.12;% parameter for re-scaling of synaptic weights
% % vars.scaleWeight = 0.12;% parameter for re-scaling of synaptic weights
% % vars.scaleWeight = 0.15;% parameter for re-scaling of synaptic weights
% % vars.scaleWeight = 0.2;% parameter for re-scaling of synaptic weights
% % vars.p_E = 0.2;% connection probability via electrical synapses
% % vars.tau_r_E = 1;% rise time for E synapses (ms)
% % vars.tau_d_E = 6;% decay time for E synapses (ms)
% % vars.tau_r_I = 0.5;% rise time for I synapses (ms)
% % vars.tau_d_I = 2;% decay time for I synapses (ms)
% vars.r_ext_EE = 12 * 500;% Rate of external input to E neurons (Hz)
% vars.r_ext_IE = 0 * 500;% Rate of external input to I neurons (Hz)
% vars.J_EE_min = 1.78;% minimum E to E synaptic weight (pF)
% vars.J_EE_max = 21.4;% maximum E to E synaptic weight (pF)
% vars.J_EE_O = 2.76;% initial E to E synaptic weight (pF)
% vars.J_EE_O1 = 2.76 * 8;% initial E to E synaptic weight (pF)
% vars.J_EI_min = 48.7;% minimum I to E synaptic weight (pF)
% vars.J_EI_max = 243;% maximum I to E synaptic weight (pF)
% vars.J_EI_O = 3 * 48.7;% initial I to E synaptic weight (pF)
% vars.J_IE_val = 1.27;% value of E to I synaptic weight (pF)
% vars.J_II_val = 16.2;% value of I to I synaptic weight (pF)
% % vars.J_GII_val = 0;% value of I to I gap junction weight (pF)
% vars.J_ext_EE = vars.J_EE_min;
% vars.J_ext_IE = vars.J_IE_val;
% vars.period_scaling = 20;% period of synaptic scaling (ms)
%
% vars.A_LTD = 0.0008;% Long-term depression (LTD) strength (pA mV^-1)
% vars.A_LTP = 0.0014;% Long-term potentiation (LTP) strength (pA mV^-2)
% vars.theta_LTD = -70;% Threshold to recruit LTD (mV)
% vars.theta_LTP = -49;% Threshold to recruit LTP (mV)
% vars.tau_u = 10;% Time constant of low-pass filtered membrane voltage (for LTD) (ms)
% vars.tau_v = 7;% Time constant of low-pass filtered membrane voltage (for LTP) (ms)
% vars.tau_x = 15;% Time constant low-pass filtered spike train (for LTP) (ms)
%
% vars.tau_y = 20;% Time constant of low-pass filtered spike train (ms)
% vars.eta = 1;% Synaptic plasticity learning rate (pA)
% vars.r_O = 3;% Target firing rate (Hz)
%
% % vars.triplet_tau_plus = 16.8;% parameter of triplet plasticity (ms)
% % vars.triplet_tau_minus = 33.7;% parameter of triplet plasticity (ms)
% % vars.triplet_tau_x = 101;% parameter of triplet plasticity (ms)
% % vars.triplet_tau_y = 125;% parameter of triplet plasticity (ms)
% % vars.triplet_A_2_plus = 7.5 * 10 ^ (-10);% parameter of triplet plasticity
% % vars.triplet_A_3_plus = 9.3 * 10 ^ (-3);% parameter of triplet plasticity
% % vars.triplet_A_2_minus = 7 * 10 ^ (-3);% parameter of triplet plasticity
% % vars.triplet_A_3_minus = 2.3 * 10 ^ (-4);% parameter of triplet plasticity
%
% vars.V_E = [];% voltages of E (mV)
% vars.V_I = [];% voltages of I (mV)
% vars.Vu_E = [];% low-pass filtered valtage of E (mV)
% vars.Vv_E = [];% low-pass filtered valtage of E (mV)
%
% vars.s_E = [];% spike train of E
% vars.sx_E = [];% low-pass filtered spike train of E
% vars.sy_E = [];% low-pass filtered spike train of E
% vars.s_I = [];% spike train of I
% vars.sy_I = [];% low-pass filtered spike train of I
%
% vars.V_T_E = [];% thereshold of E
% vars.ww_E = [];% adaptation current of E
%
% vars.s_ext_EE = [];% external stimulation to E
% vars.s_ext_IE = [];% external stimulation to I
%
% vars.ct_abs_E = [];% count of absolute refractory period for E neurons
% vars.ct_abs_I = [];% count of absolute refractory period for I neurons
%
% % vars.triplet_r_1 = [];% variable of triplet plasticity
% % vars.triplet_r_2 = [];% variable of triplet plasticity
% % vars.triplet_o_1 = [];% variable of triplet plasticity
% % vars.triplet_o_2 = [];% variable of triplet plasticity
%
% vars.connScale = [ 1, 1, 7 ];% [ orthodromic, antidromic, within area ]
% vars.connects = zeros( length( vars.N_E ), length( vars.N_E ) );% vars.connects( j, i ): from area i to area j.
% % vars.connects( 2, 1 ) = vars.connScale( 1 ) * vars.N_E( 1 );
% % vars.connects( 1, 2 ) = vars.connScale( 2 ) * vars.N_E( 2 );
% % vars.connects( 1, 1 ) = vars.connScale( 3 ) * vars.N_E( 1 );
% % vars.connects( 2, 2 ) = vars.connScale( 3 ) * vars.N_E( 2 );
% % vars.connects = vars.connects ./ sum( vars.connects, 2 );
% % vars.connects( isnan( vars.connects ) ) = 0;
% % vars.connects( vars.connects > 0 ) = 1;
% vars.connects( 2, 1 ) = vars.connScale( 1 );
% vars.connects( 3, 2 ) = vars.connScale( 1 );
% vars.connects( 1, 3 ) = vars.connScale( 2 );
% vars.connects( 1, 1 ) = vars.connScale( 3 );
% vars.connects( 2, 2 ) = vars.connScale( 3 );
% vars.connects( 3, 3 ) = vars.connScale( 3 );
%
% vars.tau_r_E = 1 * 2;% rise time for E synapses (ms)
% vars.tau_d_E = 6 * 2;% decay time for E synapses (ms)
% vars.tau_r_I = 0.5 * 10;% rise time for I synapses (ms)
% vars.tau_d_I = 2 * 10;% decay time for I synapses (ms)
% vars.nChannels = 10;
% vars = SNN_initiation0( vars );
%
%
% save( 'SNN_Simul3_0.mat', 'vars' )
% clear vars
%
%
% % dt = vars.dt;
% % t_memory = vars.t_memory;
% % t_memory_L = t_memory / dt;
% % tau_r_E = vars.tau_r_E;
% % tau_d_E = vars.tau_d_E;
% % tau_r_I = vars.tau_r_I;
% % tau_d_I = vars.tau_d_I;
% % t_domain_memory = 0 : dt : t_memory - dt;
% % F_E = ( 1 / ( tau_d_E - tau_r_E ) ) * ( exp( -t_domain_memory / tau_d_E ) - exp( -t_domain_memory / tau_r_E ) );
% % F_I = ( 1 / ( tau_d_I - tau_r_I ) ) * ( exp( -t_domain_memory / tau_d_I ) - exp( -t_domain_memory / tau_r_I ) );
% % figure; hold on; plot(F_E); plot(F_I)
%

%% SNN_Simul3_1
%
% transMatU = zeros( 30, 1500 );
% ct_u = 0;
% for u = 1 : 30
%     transMatU( u, ct_u + [ 1 : 50 ] ) = 1;
%     ct_u = ct_u + 50;
% end; clear u ct_u
%
% transMatT = zeros( 9600, 480 );
% ct_t = 0;
% for t = 1 : 480
%     transMatT( ct_t + [ 1 : 20 ], t ) = 1;
%     ct_t = ct_t + 20;
% end; clear t ct_t
%
% n_stim = 10;
% S = eye( n_stim );
%
%
% pool_E = [ 7 : 2 : 25 ];
% pool_I = [ 3 : 3 : 30 ];
%
% M = [];
% C = [];
% H = [];
% s_E_all = {};
% s_E_all_trans = {};
%
% ct_E = 0;
% for mE = pool_E
%     ct_E = ct_E + 1;
%     ct_I = 0;
%     for mI = pool_I
%         ct_I = ct_I + 1;
%
%         load( 'SNN_Simul3_0.mat' )
%
%         vars.tau_r_E = 1 * mE;% rise time for E synapses (ms)
%         vars.tau_d_E = 6 * mE * 0.3;% decay time for E synapses (ms)
%         vars.tau_r_I = 0.5 * mI;% rise time for I synapses (ms)
%         vars.tau_d_I = 2 * mI;% decay time for I synapses (ms)
%
%         vars.mode = 0.5;% the current balance between synaptic weights ( 1:= lateral, 0:= recurrent ), \in [ 0, 1 ].
%         vars.infoStims = [ 4800, 120, 40, 10, 1 ];% recording duration (ms), stimulation duration (ms), stimulation delay (ms), the number of assemblies, the number of stimulations ];
%         vars.inputs = [ 1, 2, 3 ];% in area
%         vars = SNN_stimulation( vars, 0 );
%         [ vars, records ] = SNN_simulation( vars, 'plasticity_off', 'records_on' );
%
%         s_E = records.s_E;
%         s_E_all{ ct_E, ct_I } = s_E;
%         s_E = ( transMatU * s_E ) * transMatT;
%         s_E_all_trans{ ct_E, ct_I } = s_E;
%
%         % figure
%         % imagesc( vars.ext_stim_E ); colorbar
%         % figure
%         % imagesc( records.s_E ); colorbar
%         % figure
%         % imagesc( s_E ); colorbar
%
%
%         Z1 = s_E( 1 : 10, : );
%         Z2 = s_E( 11 : 20, : );
%         Z3 = s_E( 21 : 30, : );
%
%         M1 = mean( sum( Z1, 1 ), 2 );
%         M2 = mean( sum( Z2, 1 ), 2 );
%         M3 = mean( sum( Z3, 1 ), 2 );
%         [ ~, C1, H1 ] = patternMatching( S, Z1 );
%         [ ~, C2, H2 ] = patternMatching( S, Z2 );
%         [ ~, C3, H3 ] = patternMatching( S, Z3 );
%         C1 = C1( H1 > 0 );
%         C2 = C2( H2 > 0 );
%         C3 = C3( H3 > 0 );
%         if length( C1 ) >= 0.5 * length( H1 )
%             t1C = C1( 1 : end - 1 );
%             t2C = C1( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C1 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C1 = false;
%         end
%         if length( C2 ) >= 0.5 * length( H2 )
%             t1C = C2( 1 : end - 1 );
%             t2C = C2( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C2 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C2 = false;
%         end
%         if length( C3 ) >= 0.5 * length( H2 )
%             t1C = C3( 1 : end - 1 );
%             t2C = C3( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C3 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C3 = false;
%         end
%         H1 = mean( H1, 2 );
%         H2 = mean( H2, 2 );
%         H3 = mean( H3, 2 );
%
%         M( ct_E, ct_I ) = mean( [ M1, M2, M3 ], 2 );
%         C( ct_E, ct_I ) = all( [ C1, C2, C3 ], 2 );
%         H( ct_E, ct_I ) = mean( [ H1, H2, H3 ], 2 );
%
%         disp( [ 'E: ', num2str( ct_E ), ' / I: ', num2str( ct_I ) ] )
%     end; clear mI ct_I
% end; clear mE ct_E
%
% save( 'SNN_Simul3_1.mat', 'M', 'C', 'H', 's_E_all', 's_E_all_trans', 'pool_E', 'pool_I', '-v7.3' )
%

%% Result, SNN_Simul3
%
% load( 'SNN_Simul3_0.mat' )
% load( 'SNN_Simul3_1.mat' )
%
%
% figure( 'position', [ 100, 100, 300, 220 ] )
% imagesc( vars.J_EE, [ 0, 20 ] )
% colorbar% (pF)
% set( gca, 'xtick', [], 'ytick', [] )
%
% % ct_f = 0;
% % figure
% % for ct_E = 1 : 10
% %     for ct_I = 1 : 10
% %         ct_f = ct_f + 1;
% %         subplot( 10, 10, ct_f )
% %         imagesc( s_E_all_trans{ ct_E, ct_I } )
% %     end; clear ct_I
% % end; clear ct_E
%
%
% transMatU = zeros( 30, 1500 );
% ct_u = 0;
% for u = 1 : 30
%     transMatU( u, ct_u + [ 1 : 50 ] ) = 1;
%     ct_u = ct_u + 50;
% end; clear u ct_u
%
% spkT1 = [];% first spike timing
% for ct_E = 1 : length( pool_E )
%     for ct_I = 1 : length( pool_I )
%         s_E = s_E_all{ ct_E, ct_I };
%         s_E = ( transMatU * s_E );
%         spkT1( ct_E, ct_I ) = find( s_E( 1, : ) > 25, 1 );
%         clear s_E
%     end; clear ct_I
% end; clear ct_E
% spkT1 = mean( spkT1, 2 ) / 2;% (ms)
%
%
% colors = [ 1, 1, 1; 0.8, 0.8, 0.3 ];
%
% figure( 'position', [ 100, 100, 300, 200 ] )
% imagesc( 2 * ( 0.5 + 2 ) * pool_I, 1 * pool_E, C, [ -0.5, 1.5 ] )
% colormap( gca, colors )
% colorbar( 'ticks', [ 0, 1 ], 'ticklabels', { '', '' } )
% axis xy
% hold on
% ylabel( 'Excitatory rise param. (ms)' )
% xlabel( 'Approx. inhibition dur. (ms)' )
% title( 'States of parameters' )
%
% lw = 2;
% figure( 'position', [ 100, 100, 200, 200 ] )
% plot( pool_E, spkT1, '-', 'linewidth', lw )
% xlabel( 'Excitatory rise param. (ms)' )
% ylabel( 'Actual initiation dur. (ms)' )
%
% figure( 'position', [ 100, 100, 300, 200 ] )
% imagesc( 2 * ( 0.5 + 2 ) * pool_I, spkT1, C, [ -0.5, 1.5 ] )
% colormap( gca, colors )
% colorbar( 'ticks', [ 0, 1 ], 'ticklabels', { '', '' } )
% axis xy
% hold on
% ylabel( 'Actual initiation dur. (ms)' )
% xlabel( 'Approx. inhibition dur. (ms)' )
% title( 'States of parameters' )
%

%% SNN_Simul4_0
%
% clear vars
%
% vars.dt = 0.5;% time resolution (ms)
% vars.t_memory = 250;% time memory (ms)
%
% vars.tau_E = 20;% E neuron resting membrane time constant (ms)
% vars.tau_I = 20;% I neuron resting membrane time constant (ms)
% vars.E_L_E = -70;% E neuron resting potential (mV)
% vars.E_L_I = -62;% I neuron resting potential (mV)
% vars.Delta_T_E = 2;% E neuron EIF slope factor (mV)
% vars.Delta_T_I = 2;% I neuron EIF slope factor (mV)
% vars.C = 300;% Capacitance (pF)
% vars.E_E = 0;% E reverse potential (mV)
% vars.E_I = -75;% I reverse potential (mV)
% vars.V_T = -52;% threshold potential (mV)
% vars.A_T = 10;% Post spike threshold potential increase (mV)
% vars.tau_T = 30;% Adaptive threshold time scale (ms)
% vars.V_re = -60;% Reset potential (mV)
% vars.tau_abs = 1;% Absolute refractory period (ms)
% vars.a_w = 4;% Subthreshold adaptation (nS)
% vars.b_w = 0.805;% Spike-triggered adaptation (pA)
% vars.tau_w = 150;% Spike-triggered adaptation time scale (ms)
% vars.V_ap = 20;% action potential threshold potention (mV)
% vars.V_lb = -75;% lower bound potential (mV) := vars.E_I
%
% vars.N_E = [ 500, 500, 500, 500 ];% the number of excitatory neurons
% vars.N_I = [ 125, 125, 125, 125 ];% the number of excitatory neurons
% % vars.N_E = [ 1000, 1000 ];% the number of excitatory neurons
% % vars.N_I = [ 250, 250 ];% the number of excitatory neurons
% % vars.p_C = 2.5 * 0.2;% connection probability via chemical synapses
% vars.p_C = 0.2;% connection probability via chemical synapses
% vars.scaleWeight = 0.08;% parameter for re-scaling of synaptic weights
% % vars.scaleWeight = 0.12;% parameter for re-scaling of synaptic weights
% % vars.scaleWeight = 0.15;% parameter for re-scaling of synaptic weights
% % vars.scaleWeight = 0.2;% parameter for re-scaling of synaptic weights
% % vars.p_E = 0.2;% connection probability via electrical synapses
% % vars.tau_r_E = 1;% rise time for E synapses (ms)
% % vars.tau_d_E = 6;% decay time for E synapses (ms)
% % vars.tau_r_I = 0.5;% rise time for I synapses (ms)
% % vars.tau_d_I = 2;% decay time for I synapses (ms)
% vars.r_ext_EE = 12 * 500;% Rate of external input to E neurons (Hz)
% vars.r_ext_IE = 0 * 500;% Rate of external input to I neurons (Hz)
% vars.J_EE_min = 1.78;% minimum E to E synaptic weight (pF)
% vars.J_EE_max = 21.4;% maximum E to E synaptic weight (pF)
% vars.J_EE_O = 2.76;% initial E to E synaptic weight (pF)
% vars.J_EE_O1 = 2.76 * 8;% initial E to E synaptic weight (pF)
% vars.J_EI_min = 48.7;% minimum I to E synaptic weight (pF)
% vars.J_EI_max = 243;% maximum I to E synaptic weight (pF)
% vars.J_EI_O = 3 * 48.7;% initial I to E synaptic weight (pF)
% vars.J_IE_val = 1.27;% value of E to I synaptic weight (pF)
% vars.J_II_val = 16.2;% value of I to I synaptic weight (pF)
% % vars.J_GII_val = 0;% value of I to I gap junction weight (pF)
% vars.J_ext_EE = vars.J_EE_min;
% vars.J_ext_IE = vars.J_IE_val;
% vars.period_scaling = 20;% period of synaptic scaling (ms)
%
% vars.A_LTD = 0.0008;% Long-term depression (LTD) strength (pA mV^-1)
% vars.A_LTP = 0.0014;% Long-term potentiation (LTP) strength (pA mV^-2)
% vars.theta_LTD = -70;% Threshold to recruit LTD (mV)
% vars.theta_LTP = -49;% Threshold to recruit LTP (mV)
% vars.tau_u = 10;% Time constant of low-pass filtered membrane voltage (for LTD) (ms)
% vars.tau_v = 7;% Time constant of low-pass filtered membrane voltage (for LTP) (ms)
% vars.tau_x = 15;% Time constant low-pass filtered spike train (for LTP) (ms)
%
% vars.tau_y = 20;% Time constant of low-pass filtered spike train (ms)
% vars.eta = 1;% Synaptic plasticity learning rate (pA)
% vars.r_O = 3;% Target firing rate (Hz)
%
% % vars.triplet_tau_plus = 16.8;% parameter of triplet plasticity (ms)
% % vars.triplet_tau_minus = 33.7;% parameter of triplet plasticity (ms)
% % vars.triplet_tau_x = 101;% parameter of triplet plasticity (ms)
% % vars.triplet_tau_y = 125;% parameter of triplet plasticity (ms)
% % vars.triplet_A_2_plus = 7.5 * 10 ^ (-10);% parameter of triplet plasticity
% % vars.triplet_A_3_plus = 9.3 * 10 ^ (-3);% parameter of triplet plasticity
% % vars.triplet_A_2_minus = 7 * 10 ^ (-3);% parameter of triplet plasticity
% % vars.triplet_A_3_minus = 2.3 * 10 ^ (-4);% parameter of triplet plasticity
%
% vars.V_E = [];% voltages of E (mV)
% vars.V_I = [];% voltages of I (mV)
% vars.Vu_E = [];% low-pass filtered valtage of E (mV)
% vars.Vv_E = [];% low-pass filtered valtage of E (mV)
%
% vars.s_E = [];% spike train of E
% vars.sx_E = [];% low-pass filtered spike train of E
% vars.sy_E = [];% low-pass filtered spike train of E
% vars.s_I = [];% spike train of I
% vars.sy_I = [];% low-pass filtered spike train of I
%
% vars.V_T_E = [];% thereshold of E
% vars.ww_E = [];% adaptation current of E
%
% vars.s_ext_EE = [];% external stimulation to E
% vars.s_ext_IE = [];% external stimulation to I
%
% vars.ct_abs_E = [];% count of absolute refractory period for E neurons
% vars.ct_abs_I = [];% count of absolute refractory period for I neurons
%
% % vars.triplet_r_1 = [];% variable of triplet plasticity
% % vars.triplet_r_2 = [];% variable of triplet plasticity
% % vars.triplet_o_1 = [];% variable of triplet plasticity
% % vars.triplet_o_2 = [];% variable of triplet plasticity
%
% vars.connScale = [ 1, 1, 7 ];% [ orthodromic, antidromic, within area ]
% vars.connects = zeros( length( vars.N_E ), length( vars.N_E ) );% vars.connects( j, i ): from area i to area j.
% % vars.connects( 2, 1 ) = vars.connScale( 1 ) * vars.N_E( 1 );
% % vars.connects( 1, 2 ) = vars.connScale( 2 ) * vars.N_E( 2 );
% % vars.connects( 1, 1 ) = vars.connScale( 3 ) * vars.N_E( 1 );
% % vars.connects( 2, 2 ) = vars.connScale( 3 ) * vars.N_E( 2 );
% % vars.connects = vars.connects ./ sum( vars.connects, 2 );
% % vars.connects( isnan( vars.connects ) ) = 0;
% % vars.connects( vars.connects > 0 ) = 1;
% vars.connects( 2, 1 ) = vars.connScale( 1 );
% vars.connects( 3, 2 ) = vars.connScale( 1 );
% vars.connects( 4, 3 ) = vars.connScale( 1 );
% vars.connects( 1, 4 ) = vars.connScale( 2 );
% vars.connects( 1, 1 ) = vars.connScale( 3 );
% vars.connects( 2, 2 ) = vars.connScale( 3 );
% vars.connects( 3, 3 ) = vars.connScale( 3 );
% vars.connects( 4, 4 ) = vars.connScale( 3 );
%
% vars.tau_r_E = 1 * 2;% rise time for E synapses (ms)
% vars.tau_d_E = 6 * 2;% decay time for E synapses (ms)
% vars.tau_r_I = 0.5 * 10;% rise time for I synapses (ms)
% vars.tau_d_I = 2 * 10;% decay time for I synapses (ms)
% vars.nChannels = 10;
% vars = SNN_initiation0( vars );
%
%
% save( 'SNN_Simul4_0.mat', 'vars' )
% clear vars
%
%
% % dt = vars.dt;
% % t_memory = vars.t_memory;
% % t_memory_L = t_memory / dt;
% % tau_r_E = vars.tau_r_E;
% % tau_d_E = vars.tau_d_E;
% % tau_r_I = vars.tau_r_I;
% % tau_d_I = vars.tau_d_I;
% % t_domain_memory = 0 : dt : t_memory - dt;
% % F_E = ( 1 / ( tau_d_E - tau_r_E ) ) * ( exp( -t_domain_memory / tau_d_E ) - exp( -t_domain_memory / tau_r_E ) );
% % F_I = ( 1 / ( tau_d_I - tau_r_I ) ) * ( exp( -t_domain_memory / tau_d_I ) - exp( -t_domain_memory / tau_r_I ) );
% % figure; hold on; plot(F_E); plot(F_I)
%

%% SNN_Simul4_1
%
% transMatU = zeros( 40, 2000 );
% ct_u = 0;
% for u = 1 : 40
%     transMatU( u, ct_u + [ 1 : 50 ] ) = 1;
%     ct_u = ct_u + 50;
% end; clear u ct_u
%
% transMatT = zeros( 9600, 480 );
% ct_t = 0;
% for t = 1 : 480
%     transMatT( ct_t + [ 1 : 20 ], t ) = 1;
%     ct_t = ct_t + 20;
% end; clear t ct_t
%
% n_stim = 10;
% S = eye( n_stim );
%
%
% pool_E = [ 7 : 2 : 25 ];
% pool_I = [ 3 : 3 : 30 ];
%
% M = [];
% C = [];
% H = [];
% s_E_all = {};
% s_E_all_trans = {};
%
% ct_E = 0;
% for mE = pool_E
%     ct_E = ct_E + 1;
%     ct_I = 0;
%     for mI = pool_I
%         ct_I = ct_I + 1;
%
%         load( 'SNN_Simul4_0.mat' )
%
%         vars.tau_r_E = 1 * mE;% rise time for E synapses (ms)
%         vars.tau_d_E = 6 * mE * 0.3;% decay time for E synapses (ms)
%         vars.tau_r_I = 0.5 * mI;% rise time for I synapses (ms)
%         vars.tau_d_I = 2 * mI;% decay time for I synapses (ms)
%
%         vars.mode = 0.5;% the current balance between synaptic weights ( 1:= lateral, 0:= recurrent ), \in [ 0, 1 ].
%         vars.infoStims = [ 4800, 120, 40, 10, 1 ];% recording duration (ms), stimulation duration (ms), stimulation delay (ms), the number of assemblies, the number of stimulations ];
%         vars.inputs = [ 1, 2, 3, 4 ];% in area
%         vars = SNN_stimulation( vars, 0 );
%         [ vars, records ] = SNN_simulation( vars, 'plasticity_off', 'records_on' );
%
%         s_E = records.s_E;
%         s_E_all{ ct_E, ct_I } = s_E;
%         s_E = ( transMatU * s_E ) * transMatT;
%         s_E_all_trans{ ct_E, ct_I } = s_E;
%
%         % figure
%         % imagesc( vars.ext_stim_E ); colorbar
%         % figure
%         % imagesc( records.s_E ); colorbar
%         % figure
%         % imagesc( s_E ); colorbar
%
%
%         Z1 = s_E( 1 : 10, : );
%         Z2 = s_E( 11 : 20, : );
%         Z3 = s_E( 21 : 30, : );
%         Z4 = s_E( 31 : 40, : );
%
%         M1 = mean( sum( Z1, 1 ), 2 );
%         M2 = mean( sum( Z2, 1 ), 2 );
%         M3 = mean( sum( Z3, 1 ), 2 );
%         M4 = mean( sum( Z4, 1 ), 2 );
%         [ ~, C1, H1 ] = patternMatching( S, Z1 );
%         [ ~, C2, H2 ] = patternMatching( S, Z2 );
%         [ ~, C3, H3 ] = patternMatching( S, Z3 );
%         [ ~, C4, H4 ] = patternMatching( S, Z4 );
%         C1 = C1( H1 > 0 );
%         C2 = C2( H2 > 0 );
%         C3 = C3( H3 > 0 );
%         C4 = C4( H4 > 0 );
%         if length( C1 ) >= 0.5 * length( H1 )
%             t1C = C1( 1 : end - 1 );
%             t2C = C1( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C1 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C1 = false;
%         end
%         if length( C2 ) >= 0.5 * length( H2 )
%             t1C = C2( 1 : end - 1 );
%             t2C = C2( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C2 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C2 = false;
%         end
%         if length( C3 ) >= 0.5 * length( H2 )
%             t1C = C3( 1 : end - 1 );
%             t2C = C3( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C3 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C3 = false;
%         end
%         if length( C4 ) >= 0.5 * length( H2 )
%             t1C = C4( 1 : end - 1 );
%             t2C = C4( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C4 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C4 = false;
%         end
%         H1 = mean( H1, 2 );
%         H2 = mean( H2, 2 );
%         H3 = mean( H3, 2 );
%         H4 = mean( H4, 2 );
%
%         M( ct_E, ct_I ) = mean( [ M1, M2, M3, M4 ], 2 );
%         C( ct_E, ct_I ) = all( [ C1, C2, C3, C4 ], 2 );
%         H( ct_E, ct_I ) = mean( [ H1, H2, H3, H4 ], 2 );
%
%         disp( [ 'E: ', num2str( ct_E ), ' / I: ', num2str( ct_I ) ] )
%     end; clear mI ct_I
% end; clear mE ct_E
%
% save( 'SNN_Simul4_1.mat', 'M', 'C', 'H', 's_E_all', 's_E_all_trans', 'pool_E', 'pool_I', '-v7.3' )
%

%% Result, SNN_Simul4
%
% load( 'SNN_Simul4_0.mat' )
% load( 'SNN_Simul4_1.mat' )
%
%
% figure( 'position', [ 100, 100, 300, 220 ] )
% imagesc( vars.J_EE, [ 0, 20 ] )
% colorbar% (pF)
% set( gca, 'xtick', [], 'ytick', [] )
%
% % ct_f = 0;
% % figure
% % for ct_E = 1 : 10
% %     for ct_I = 1 : 10
% %         ct_f = ct_f + 1;
% %         subplot( 10, 10, ct_f )
% %         imagesc( s_E_all_trans{ ct_E, ct_I } )
% %     end; clear ct_I
% % end; clear ct_E
%
%
% transMatU = zeros( 40, 2000 );
% ct_u = 0;
% for u = 1 : 40
%     transMatU( u, ct_u + [ 1 : 50 ] ) = 1;
%     ct_u = ct_u + 50;
% end; clear u ct_u
%
% spkT1 = [];% first spike timing
% for ct_E = 1 : length( pool_E )
%     for ct_I = 1 : length( pool_I )
%         s_E = s_E_all{ ct_E, ct_I };
%         s_E = ( transMatU * s_E );
%         spkT1( ct_E, ct_I ) = find( s_E( 1, : ) > 25, 1 );
%         clear s_E
%     end; clear ct_I
% end; clear ct_E
% spkT1 = mean( spkT1, 2 ) / 2;% (ms)
%
%
% colors = [ 1, 1, 1; 0.8, 0.8, 0.3 ];
%
% figure( 'position', [ 100, 100, 300, 200 ] )
% imagesc( 2 * ( 0.5 + 2 ) * pool_I, 1 * pool_E, C, [ -0.5, 1.5 ] )
% colormap( gca, colors )
% colorbar( 'ticks', [ 0, 1 ], 'ticklabels', { '', '' } )
% axis xy
% hold on
% ylabel( 'Excitatory rise param. (ms)' )
% xlabel( 'Approx. inhibition dur. (ms)' )
% title( 'States of parameters' )
%
% lw = 2;
% figure( 'position', [ 100, 100, 200, 200 ] )
% plot( pool_E, spkT1, '-', 'linewidth', lw )
% xlabel( 'Excitatory rise param. (ms)' )
% ylabel( 'Actual initiation dur. (ms)' )
%
% figure( 'position', [ 100, 100, 300, 200 ] )
% imagesc( 2 * ( 0.5 + 2 ) * pool_I, spkT1, C, [ -0.5, 1.5 ] )
% colormap( gca, colors )
% colorbar( 'ticks', [ 0, 1 ], 'ticklabels', { '', '' } )
% axis xy
% hold on
% ylabel( 'Actual initiation dur. (ms)' )
% xlabel( 'Approx. inhibition dur. (ms)' )
% title( 'States of parameters' )
%

%% SNN_Simul5_0
%
% clear vars
%
% vars.dt = 0.5;% time resolution (ms)
% vars.t_memory = 250;% time memory (ms)
%
% vars.tau_E = 20;% E neuron resting membrane time constant (ms)
% vars.tau_I = 20;% I neuron resting membrane time constant (ms)
% vars.E_L_E = -70;% E neuron resting potential (mV)
% vars.E_L_I = -62;% I neuron resting potential (mV)
% vars.Delta_T_E = 2;% E neuron EIF slope factor (mV)
% vars.Delta_T_I = 2;% I neuron EIF slope factor (mV)
% vars.C = 300;% Capacitance (pF)
% vars.E_E = 0;% E reverse potential (mV)
% vars.E_I = -75;% I reverse potential (mV)
% vars.V_T = -52;% threshold potential (mV)
% vars.A_T = 10;% Post spike threshold potential increase (mV)
% vars.tau_T = 30;% Adaptive threshold time scale (ms)
% vars.V_re = -60;% Reset potential (mV)
% vars.tau_abs = 1;% Absolute refractory period (ms)
% vars.a_w = 4;% Subthreshold adaptation (nS)
% vars.b_w = 0.805;% Spike-triggered adaptation (pA)
% vars.tau_w = 150;% Spike-triggered adaptation time scale (ms)
% vars.V_ap = 20;% action potential threshold potention (mV)
% vars.V_lb = -75;% lower bound potential (mV) := vars.E_I
%
% vars.N_E = [ 500, 500, 500, 500, 500 ];% the number of excitatory neurons
% vars.N_I = [ 125, 125, 125, 125, 125 ];% the number of excitatory neurons
% % vars.N_E = [ 1000, 1000 ];% the number of excitatory neurons
% % vars.N_I = [ 250, 250 ];% the number of excitatory neurons
% % vars.p_C = 2.5 * 0.2;% connection probability via chemical synapses
% vars.p_C = 0.2;% connection probability via chemical synapses
% vars.scaleWeight = 0.06;% parameter for re-scaling of synaptic weights
% % vars.scaleWeight = 0.12;% parameter for re-scaling of synaptic weights
% % vars.scaleWeight = 0.15;% parameter for re-scaling of synaptic weights
% % vars.scaleWeight = 0.2;% parameter for re-scaling of synaptic weights
% % vars.p_E = 0.2;% connection probability via electrical synapses
% % vars.tau_r_E = 1;% rise time for E synapses (ms)
% % vars.tau_d_E = 6;% decay time for E synapses (ms)
% % vars.tau_r_I = 0.5;% rise time for I synapses (ms)
% % vars.tau_d_I = 2;% decay time for I synapses (ms)
% vars.r_ext_EE = 12 * 500;% Rate of external input to E neurons (Hz)
% vars.r_ext_IE = 0 * 500;% Rate of external input to I neurons (Hz)
% vars.J_EE_min = 1.78;% minimum E to E synaptic weight (pF)
% vars.J_EE_max = 21.4;% maximum E to E synaptic weight (pF)
% vars.J_EE_O = 2.76;% initial E to E synaptic weight (pF)
% vars.J_EE_O1 = 2.76 * 8;% initial E to E synaptic weight (pF)
% vars.J_EI_min = 48.7;% minimum I to E synaptic weight (pF)
% vars.J_EI_max = 243;% maximum I to E synaptic weight (pF)
% vars.J_EI_O = 3 * 48.7;% initial I to E synaptic weight (pF)
% vars.J_IE_val = 1.27;% value of E to I synaptic weight (pF)
% vars.J_II_val = 16.2;% value of I to I synaptic weight (pF)
% % vars.J_GII_val = 0;% value of I to I gap junction weight (pF)
% vars.J_ext_EE = vars.J_EE_min;
% vars.J_ext_IE = vars.J_IE_val;
% vars.period_scaling = 20;% period of synaptic scaling (ms)
%
% vars.A_LTD = 0.0008;% Long-term depression (LTD) strength (pA mV^-1)
% vars.A_LTP = 0.0014;% Long-term potentiation (LTP) strength (pA mV^-2)
% vars.theta_LTD = -70;% Threshold to recruit LTD (mV)
% vars.theta_LTP = -49;% Threshold to recruit LTP (mV)
% vars.tau_u = 10;% Time constant of low-pass filtered membrane voltage (for LTD) (ms)
% vars.tau_v = 7;% Time constant of low-pass filtered membrane voltage (for LTP) (ms)
% vars.tau_x = 15;% Time constant low-pass filtered spike train (for LTP) (ms)
%
% vars.tau_y = 20;% Time constant of low-pass filtered spike train (ms)
% vars.eta = 1;% Synaptic plasticity learning rate (pA)
% vars.r_O = 3;% Target firing rate (Hz)
%
% % vars.triplet_tau_plus = 16.8;% parameter of triplet plasticity (ms)
% % vars.triplet_tau_minus = 33.7;% parameter of triplet plasticity (ms)
% % vars.triplet_tau_x = 101;% parameter of triplet plasticity (ms)
% % vars.triplet_tau_y = 125;% parameter of triplet plasticity (ms)
% % vars.triplet_A_2_plus = 7.5 * 10 ^ (-10);% parameter of triplet plasticity
% % vars.triplet_A_3_plus = 9.3 * 10 ^ (-3);% parameter of triplet plasticity
% % vars.triplet_A_2_minus = 7 * 10 ^ (-3);% parameter of triplet plasticity
% % vars.triplet_A_3_minus = 2.3 * 10 ^ (-4);% parameter of triplet plasticity
%
% vars.V_E = [];% voltages of E (mV)
% vars.V_I = [];% voltages of I (mV)
% vars.Vu_E = [];% low-pass filtered valtage of E (mV)
% vars.Vv_E = [];% low-pass filtered valtage of E (mV)
%
% vars.s_E = [];% spike train of E
% vars.sx_E = [];% low-pass filtered spike train of E
% vars.sy_E = [];% low-pass filtered spike train of E
% vars.s_I = [];% spike train of I
% vars.sy_I = [];% low-pass filtered spike train of I
%
% vars.V_T_E = [];% thereshold of E
% vars.ww_E = [];% adaptation current of E
%
% vars.s_ext_EE = [];% external stimulation to E
% vars.s_ext_IE = [];% external stimulation to I
%
% vars.ct_abs_E = [];% count of absolute refractory period for E neurons
% vars.ct_abs_I = [];% count of absolute refractory period for I neurons
%
% % vars.triplet_r_1 = [];% variable of triplet plasticity
% % vars.triplet_r_2 = [];% variable of triplet plasticity
% % vars.triplet_o_1 = [];% variable of triplet plasticity
% % vars.triplet_o_2 = [];% variable of triplet plasticity
%
% vars.connScale = [ 1, 1, 7 ];% [ orthodromic, antidromic, within area ]
% vars.connects = zeros( length( vars.N_E ), length( vars.N_E ) );% vars.connects( j, i ): from area i to area j.
% % vars.connects( 2, 1 ) = vars.connScale( 1 ) * vars.N_E( 1 );
% % vars.connects( 1, 2 ) = vars.connScale( 2 ) * vars.N_E( 2 );
% % vars.connects( 1, 1 ) = vars.connScale( 3 ) * vars.N_E( 1 );
% % vars.connects( 2, 2 ) = vars.connScale( 3 ) * vars.N_E( 2 );
% % vars.connects = vars.connects ./ sum( vars.connects, 2 );
% % vars.connects( isnan( vars.connects ) ) = 0;
% % vars.connects( vars.connects > 0 ) = 1;
% vars.connects( 2, 1 ) = vars.connScale( 1 );
% vars.connects( 3, 2 ) = vars.connScale( 1 );
% vars.connects( 4, 3 ) = vars.connScale( 1 );
% vars.connects( 5, 4 ) = vars.connScale( 1 );
% vars.connects( 1, 5 ) = vars.connScale( 2 );
% vars.connects( 1, 1 ) = vars.connScale( 3 );
% vars.connects( 2, 2 ) = vars.connScale( 3 );
% vars.connects( 3, 3 ) = vars.connScale( 3 );
% vars.connects( 4, 4 ) = vars.connScale( 3 );
% vars.connects( 5, 5 ) = vars.connScale( 3 );
%
% vars.tau_r_E = 1 * 2;% rise time for E synapses (ms)
% vars.tau_d_E = 6 * 2;% decay time for E synapses (ms)
% vars.tau_r_I = 0.5 * 10;% rise time for I synapses (ms)
% vars.tau_d_I = 2 * 10;% decay time for I synapses (ms)
% vars.nChannels = 10;
% vars = SNN_initiation0( vars );
%
%
% save( 'SNN_Simul5_0.mat', 'vars' )
% clear vars
%
%
% % dt = vars.dt;
% % t_memory = vars.t_memory;
% % t_memory_L = t_memory / dt;
% % tau_r_E = vars.tau_r_E;
% % tau_d_E = vars.tau_d_E;
% % tau_r_I = vars.tau_r_I;
% % tau_d_I = vars.tau_d_I;
% % t_domain_memory = 0 : dt : t_memory - dt;
% % F_E = ( 1 / ( tau_d_E - tau_r_E ) ) * ( exp( -t_domain_memory / tau_d_E ) - exp( -t_domain_memory / tau_r_E ) );
% % F_I = ( 1 / ( tau_d_I - tau_r_I ) ) * ( exp( -t_domain_memory / tau_d_I ) - exp( -t_domain_memory / tau_r_I ) );
% % figure; hold on; plot(F_E); plot(F_I)
%

%% SNN_Simul5_1
%
% transMatU = zeros( 50, 2500 );
% ct_u = 0;
% for u = 1 : 50
%     transMatU( u, ct_u + [ 1 : 50 ] ) = 1;
%     ct_u = ct_u + 50;
% end; clear u ct_u
%
% transMatT = zeros( 9600, 480 );
% ct_t = 0;
% for t = 1 : 480
%     transMatT( ct_t + [ 1 : 20 ], t ) = 1;
%     ct_t = ct_t + 20;
% end; clear t ct_t
%
% n_stim = 10;
% S = eye( n_stim );
%
%
% pool_E = [ 7 : 2 : 25 ];
% pool_I = [ 3 : 3 : 30 ];
%
% M = [];
% C = [];
% H = [];
% s_E_all = {};
% s_E_all_trans = {};
%
% ct_E = 0;
% for mE = pool_E
%     ct_E = ct_E + 1;
%     ct_I = 0;
%     for mI = pool_I
%         ct_I = ct_I + 1;
%
%         load( 'SNN_Simul5_0.mat' )
%
%         vars.tau_r_E = 1 * mE;% rise time for E synapses (ms)
%         vars.tau_d_E = 6 * mE * 0.3;% decay time for E synapses (ms)
%         vars.tau_r_I = 0.5 * mI;% rise time for I synapses (ms)
%         vars.tau_d_I = 2 * mI;% decay time for I synapses (ms)
%
%         vars.mode = 0.5;% the current balance between synaptic weights ( 1:= lateral, 0:= recurrent ), \in [ 0, 1 ].
%         vars.infoStims = [ 4800, 120, 40, 10, 1 ];% recording duration (ms), stimulation duration (ms), stimulation delay (ms), the number of assemblies, the number of stimulations ];
%         vars.inputs = [ 1, 2, 3, 4, 5 ];% in area
%         vars = SNN_stimulation( vars, 0 );
%         [ vars, records ] = SNN_simulation_GPU( vars, 'plasticity_off', 'records_on' );
%
%         s_E = records.s_E;
%         s_E_all{ ct_E, ct_I } = s_E;
%         s_E = ( transMatU * s_E ) * transMatT;
%         s_E_all_trans{ ct_E, ct_I } = s_E;
%
%         % figure
%         % imagesc( vars.ext_stim_E ); colorbar
%         % figure
%         % imagesc( records.s_E ); colorbar
%         % figure
%         % imagesc( s_E ); colorbar
%
%
%         Z1 = s_E( 1 : 10, : );
%         Z2 = s_E( 11 : 20, : );
%         Z3 = s_E( 21 : 30, : );
%         Z4 = s_E( 31 : 40, : );
%         Z5 = s_E( 41 : 50, : );
%
%         M1 = mean( sum( Z1, 1 ), 2 );
%         M2 = mean( sum( Z2, 1 ), 2 );
%         M3 = mean( sum( Z3, 1 ), 2 );
%         M4 = mean( sum( Z4, 1 ), 2 );
%         M5 = mean( sum( Z5, 1 ), 2 );
%         [ ~, C1, H1 ] = patternMatching( S, Z1 );
%         [ ~, C2, H2 ] = patternMatching( S, Z2 );
%         [ ~, C3, H3 ] = patternMatching( S, Z3 );
%         [ ~, C4, H4 ] = patternMatching( S, Z4 );
%         [ ~, C5, H5 ] = patternMatching( S, Z5 );
%         C1 = C1( H1 > 0 );
%         C2 = C2( H2 > 0 );
%         C3 = C3( H3 > 0 );
%         C4 = C4( H4 > 0 );
%         C5 = C5( H5 > 0 );
%         if length( C1 ) >= 0.5 * length( H1 )
%             t1C = C1( 1 : end - 1 );
%             t2C = C1( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C1 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C1 = false;
%         end
%         if length( C2 ) >= 0.5 * length( H2 )
%             t1C = C2( 1 : end - 1 );
%             t2C = C2( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C2 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C2 = false;
%         end
%         if length( C3 ) >= 0.5 * length( H2 )
%             t1C = C3( 1 : end - 1 );
%             t2C = C3( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C3 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C3 = false;
%         end
%         if length( C4 ) >= 0.5 * length( H2 )
%             t1C = C4( 1 : end - 1 );
%             t2C = C4( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C4 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C4 = false;
%         end
%         if length( C5 ) >= 0.5 * length( H2 )
%             t1C = C5( 1 : end - 1 );
%             t2C = C5( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C5 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C5 = false;
%         end
%         H1 = mean( H1, 2 );
%         H2 = mean( H2, 2 );
%         H3 = mean( H3, 2 );
%         H4 = mean( H4, 2 );
%         H5 = mean( H5, 2 );
%
%         M( ct_E, ct_I ) = mean( [ M1, M2, M3, M4, M5 ], 2 );
%         C( ct_E, ct_I ) = all( [ C1, C2, C3, C4, C5 ], 2 );
%         H( ct_E, ct_I ) = mean( [ H1, H2, H3, H4, H5 ], 2 );
%
%         disp( [ 'E: ', num2str( ct_E ), ' / I: ', num2str( ct_I ) ] )
%     end; clear mI ct_I
% end; clear mE ct_E
%
% save( 'SNN_Simul5_1.mat', 'M', 'C', 'H', 's_E_all', 's_E_all_trans', 'pool_E', 'pool_I', '-v7.3' )
%

%% Result, SNN_Simul5
%
% load( 'SNN_Simul5_0.mat' )
% load( 'SNN_Simul5_1.mat' )
%
%
% figure( 'position', [ 100, 100, 300, 220 ] )
% imagesc( vars.J_EE, [ 0, 20 ] )
% colorbar% (pF)
% set( gca, 'xtick', [], 'ytick', [] )
%
% % ct_f = 0;
% % figure
% % for ct_E = 1 : 10
% %     for ct_I = 1 : 10
% %         ct_f = ct_f + 1;
% %         subplot( 10, 10, ct_f )
% %         imagesc( s_E_all_trans{ ct_E, ct_I } )
% %     end; clear ct_I
% % end; clear ct_E
%
%
% transMatU = zeros( 50, 2500 );
% ct_u = 0;
% for u = 1 : 50
%     transMatU( u, ct_u + [ 1 : 50 ] ) = 1;
%     ct_u = ct_u + 50;
% end; clear u ct_u
%
% spkT1 = [];% first spike timing
% for ct_E = 1 : length( pool_E )
%     for ct_I = 1 : length( pool_I )
%         s_E = s_E_all{ ct_E, ct_I };
%         s_E = ( transMatU * s_E );
%         spkT1( ct_E, ct_I ) = find( s_E( 1, : ) > 25, 1 );
%         clear s_E
%     end; clear ct_I
% end; clear ct_E
% spkT1 = mean( spkT1, 2 ) / 2;% (ms)
%
%
% colors = [ 1, 1, 1; 0.8, 0.8, 0.3 ];
%
% figure( 'position', [ 100, 100, 300, 200 ] )
% imagesc( 2 * ( 0.5 + 2 ) * pool_I, 1 * pool_E, C, [ -0.5, 1.5 ] )
% colormap( gca, colors )
% colorbar( 'ticks', [ 0, 1 ], 'ticklabels', { '', '' } )
% axis xy
% hold on
% ylabel( 'Excitatory rise param. (ms)' )
% xlabel( 'Approx. inhibition dur. (ms)' )
% title( 'States of parameters' )
%
% lw = 2;
% figure( 'position', [ 100, 100, 200, 200 ] )
% plot( pool_E, spkT1, '-', 'linewidth', lw )
% xlabel( 'Excitatory rise param. (ms)' )
% ylabel( 'Actual initiation dur. (ms)' )
%
% figure( 'position', [ 100, 100, 300, 200 ] )
% imagesc( 2 * ( 0.5 + 2 ) * pool_I, spkT1, C, [ -0.5, 1.5 ] )
% colormap( gca, colors )
% colorbar( 'ticks', [ 0, 1 ], 'ticklabels', { '', '' } )
% axis xy
% hold on
% ylabel( 'Actual initiation dur. (ms)' )
% xlabel( 'Approx. inhibition dur. (ms)' )
% title( 'States of parameters' )
%

%% SNN_Simul6_0
%
% clear vars
%
% vars.dt = 0.5;% time resolution (ms)
% vars.t_memory = 250;% time memory (ms)
%
% vars.tau_E = 20;% E neuron resting membrane time constant (ms)
% vars.tau_I = 20;% I neuron resting membrane time constant (ms)
% vars.E_L_E = -70;% E neuron resting potential (mV)
% vars.E_L_I = -62;% I neuron resting potential (mV)
% vars.Delta_T_E = 2;% E neuron EIF slope factor (mV)
% vars.Delta_T_I = 2;% I neuron EIF slope factor (mV)
% vars.C = 300;% Capacitance (pF)
% vars.E_E = 0;% E reverse potential (mV)
% vars.E_I = -75;% I reverse potential (mV)
% vars.V_T = -52;% threshold potential (mV)
% vars.A_T = 10;% Post spike threshold potential increase (mV)
% vars.tau_T = 30;% Adaptive threshold time scale (ms)
% vars.V_re = -60;% Reset potential (mV)
% vars.tau_abs = 1;% Absolute refractory period (ms)
% vars.a_w = 4;% Subthreshold adaptation (nS)
% vars.b_w = 0.805;% Spike-triggered adaptation (pA)
% vars.tau_w = 150;% Spike-triggered adaptation time scale (ms)
% vars.V_ap = 20;% action potential threshold potention (mV)
% vars.V_lb = -75;% lower bound potential (mV) := vars.E_I
%
% vars.N_E = [ 500, 500, 500, 500, 500, 500 ];% the number of excitatory neurons
% vars.N_I = [ 125, 125, 125, 125, 125, 125 ];% the number of excitatory neurons
% % vars.N_E = [ 1000, 1000 ];% the number of excitatory neurons
% % vars.N_I = [ 250, 250 ];% the number of excitatory neurons
% % vars.p_C = 2.5 * 0.2;% connection probability via chemical synapses
% vars.p_C = 0.2;% connection probability via chemical synapses
% vars.scaleWeight = 0.05;% parameter for re-scaling of synaptic weights
% % vars.scaleWeight = 0.12;% parameter for re-scaling of synaptic weights
% % vars.scaleWeight = 0.15;% parameter for re-scaling of synaptic weights
% % vars.scaleWeight = 0.2;% parameter for re-scaling of synaptic weights
% % vars.p_E = 0.2;% connection probability via electrical synapses
% % vars.tau_r_E = 1;% rise time for E synapses (ms)
% % vars.tau_d_E = 6;% decay time for E synapses (ms)
% % vars.tau_r_I = 0.5;% rise time for I synapses (ms)
% % vars.tau_d_I = 2;% decay time for I synapses (ms)
% vars.r_ext_EE = 12 * 500;% Rate of external input to E neurons (Hz)
% vars.r_ext_IE = 0 * 500;% Rate of external input to I neurons (Hz)
% vars.J_EE_min = 1.78;% minimum E to E synaptic weight (pF)
% vars.J_EE_max = 21.4;% maximum E to E synaptic weight (pF)
% vars.J_EE_O = 2.76;% initial E to E synaptic weight (pF)
% vars.J_EE_O1 = 2.76 * 8;% initial E to E synaptic weight (pF)
% vars.J_EI_min = 48.7;% minimum I to E synaptic weight (pF)
% vars.J_EI_max = 243;% maximum I to E synaptic weight (pF)
% vars.J_EI_O = 3 * 48.7;% initial I to E synaptic weight (pF)
% vars.J_IE_val = 1.27;% value of E to I synaptic weight (pF)
% vars.J_II_val = 16.2;% value of I to I synaptic weight (pF)
% % vars.J_GII_val = 0;% value of I to I gap junction weight (pF)
% vars.J_ext_EE = vars.J_EE_min;
% vars.J_ext_IE = vars.J_IE_val;
% vars.period_scaling = 20;% period of synaptic scaling (ms)
%
% vars.A_LTD = 0.0008;% Long-term depression (LTD) strength (pA mV^-1)
% vars.A_LTP = 0.0014;% Long-term potentiation (LTP) strength (pA mV^-2)
% vars.theta_LTD = -70;% Threshold to recruit LTD (mV)
% vars.theta_LTP = -49;% Threshold to recruit LTP (mV)
% vars.tau_u = 10;% Time constant of low-pass filtered membrane voltage (for LTD) (ms)
% vars.tau_v = 7;% Time constant of low-pass filtered membrane voltage (for LTP) (ms)
% vars.tau_x = 15;% Time constant low-pass filtered spike train (for LTP) (ms)
%
% vars.tau_y = 20;% Time constant of low-pass filtered spike train (ms)
% vars.eta = 1;% Synaptic plasticity learning rate (pA)
% vars.r_O = 3;% Target firing rate (Hz)
%
% % vars.triplet_tau_plus = 16.8;% parameter of triplet plasticity (ms)
% % vars.triplet_tau_minus = 33.7;% parameter of triplet plasticity (ms)
% % vars.triplet_tau_x = 101;% parameter of triplet plasticity (ms)
% % vars.triplet_tau_y = 125;% parameter of triplet plasticity (ms)
% % vars.triplet_A_2_plus = 7.5 * 10 ^ (-10);% parameter of triplet plasticity
% % vars.triplet_A_3_plus = 9.3 * 10 ^ (-3);% parameter of triplet plasticity
% % vars.triplet_A_2_minus = 7 * 10 ^ (-3);% parameter of triplet plasticity
% % vars.triplet_A_3_minus = 2.3 * 10 ^ (-4);% parameter of triplet plasticity
%
% vars.V_E = [];% voltages of E (mV)
% vars.V_I = [];% voltages of I (mV)
% vars.Vu_E = [];% low-pass filtered valtage of E (mV)
% vars.Vv_E = [];% low-pass filtered valtage of E (mV)
%
% vars.s_E = [];% spike train of E
% vars.sx_E = [];% low-pass filtered spike train of E
% vars.sy_E = [];% low-pass filtered spike train of E
% vars.s_I = [];% spike train of I
% vars.sy_I = [];% low-pass filtered spike train of I
%
% vars.V_T_E = [];% thereshold of E
% vars.ww_E = [];% adaptation current of E
%
% vars.s_ext_EE = [];% external stimulation to E
% vars.s_ext_IE = [];% external stimulation to I
%
% vars.ct_abs_E = [];% count of absolute refractory period for E neurons
% vars.ct_abs_I = [];% count of absolute refractory period for I neurons
%
% % vars.triplet_r_1 = [];% variable of triplet plasticity
% % vars.triplet_r_2 = [];% variable of triplet plasticity
% % vars.triplet_o_1 = [];% variable of triplet plasticity
% % vars.triplet_o_2 = [];% variable of triplet plasticity
%
% vars.connScale = [ 1, 1, 7 ];% [ orthodromic, antidromic, within area ]
% vars.connects = zeros( length( vars.N_E ), length( vars.N_E ) );% vars.connects( j, i ): from area i to area j.
% % vars.connects( 2, 1 ) = vars.connScale( 1 ) * vars.N_E( 1 );
% % vars.connects( 1, 2 ) = vars.connScale( 2 ) * vars.N_E( 2 );
% % vars.connects( 1, 1 ) = vars.connScale( 3 ) * vars.N_E( 1 );
% % vars.connects( 2, 2 ) = vars.connScale( 3 ) * vars.N_E( 2 );
% % vars.connects = vars.connects ./ sum( vars.connects, 2 );
% % vars.connects( isnan( vars.connects ) ) = 0;
% % vars.connects( vars.connects > 0 ) = 1;
% vars.connects( 2, 1 ) = vars.connScale( 1 );
% vars.connects( 3, 2 ) = vars.connScale( 1 );
% vars.connects( 4, 3 ) = vars.connScale( 1 );
% vars.connects( 5, 4 ) = vars.connScale( 1 );
% vars.connects( 6, 5 ) = vars.connScale( 1 );
% vars.connects( 1, 6 ) = vars.connScale( 2 );
% vars.connects( 1, 1 ) = vars.connScale( 3 );
% vars.connects( 2, 2 ) = vars.connScale( 3 );
% vars.connects( 3, 3 ) = vars.connScale( 3 );
% vars.connects( 4, 4 ) = vars.connScale( 3 );
% vars.connects( 5, 5 ) = vars.connScale( 3 );
% vars.connects( 6, 6 ) = vars.connScale( 3 );
%
% vars.tau_r_E = 1 * 2;% rise time for E synapses (ms)
% vars.tau_d_E = 6 * 2;% decay time for E synapses (ms)
% vars.tau_r_I = 0.5 * 10;% rise time for I synapses (ms)
% vars.tau_d_I = 2 * 10;% decay time for I synapses (ms)
% vars.nChannels = 10;
% vars = SNN_initiation0( vars );
%
%
% save( 'SNN_Simul6_0.mat', 'vars' )
% clear vars
%
%
% % dt = vars.dt;
% % t_memory = vars.t_memory;
% % t_memory_L = t_memory / dt;
% % tau_r_E = vars.tau_r_E;
% % tau_d_E = vars.tau_d_E;
% % tau_r_I = vars.tau_r_I;
% % tau_d_I = vars.tau_d_I;
% % t_domain_memory = 0 : dt : t_memory - dt;
% % F_E = ( 1 / ( tau_d_E - tau_r_E ) ) * ( exp( -t_domain_memory / tau_d_E ) - exp( -t_domain_memory / tau_r_E ) );
% % F_I = ( 1 / ( tau_d_I - tau_r_I ) ) * ( exp( -t_domain_memory / tau_d_I ) - exp( -t_domain_memory / tau_r_I ) );
% % figure; hold on; plot(F_E); plot(F_I)
%

%% SNN_Simul6_1
%
% transMatU = zeros( 60, 3000 );
% ct_u = 0;
% for u = 1 : 60
%     transMatU( u, ct_u + [ 1 : 50 ] ) = 1;
%     ct_u = ct_u + 50;
% end; clear u ct_u
%
% transMatT = zeros( 9600, 480 );
% ct_t = 0;
% for t = 1 : 480
%     transMatT( ct_t + [ 1 : 20 ], t ) = 1;
%     ct_t = ct_t + 20;
% end; clear t ct_t
%
% n_stim = 10;
% S = eye( n_stim );
%
%
% pool_E = [ 7 : 2 : 25 ];
% pool_I = [ 3 : 3 : 30 ];
%
% M = [];
% C = [];
% H = [];
% s_E_all = {};
% s_E_all_trans = {};
%
% ct_E = 0;
% for mE = pool_E
%     ct_E = ct_E + 1;
%     ct_I = 0;
%     for mI = pool_I
%         ct_I = ct_I + 1;
%
%         load( 'SNN_Simul6_0.mat' )
%
%         vars.tau_r_E = 1 * mE;% rise time for E synapses (ms)
%         vars.tau_d_E = 6 * mE * 0.3;% decay time for E synapses (ms)
%         vars.tau_r_I = 0.5 * mI;% rise time for I synapses (ms)
%         vars.tau_d_I = 2 * mI;% decay time for I synapses (ms)
%
%         vars.mode = 0.5;% the current balance between synaptic weights ( 1:= lateral, 0:= recurrent ), \in [ 0, 1 ].
%         vars.infoStims = [ 4800, 120, 40, 10, 1 ];% recording duration (ms), stimulation duration (ms), stimulation delay (ms), the number of assemblies, the number of stimulations ];
%         vars.inputs = [ 1, 2, 3, 4, 5, 6 ];% in area
%         vars = SNN_stimulation( vars, 0 );
%         [ vars, records ] = SNN_simulation_GPU( vars, 'plasticity_off', 'records_on' );
%
%         s_E = records.s_E;
%         s_E_all{ ct_E, ct_I } = s_E;
%         s_E = ( transMatU * s_E ) * transMatT;
%         s_E_all_trans{ ct_E, ct_I } = s_E;
%
%         % figure
%         % imagesc( vars.ext_stim_E ); colorbar
%         % figure
%         % imagesc( records.s_E ); colorbar
%         % figure
%         % imagesc( s_E ); colorbar
%
%
%         Z1 = s_E( 1 : 10, : );
%         Z2 = s_E( 11 : 20, : );
%         Z3 = s_E( 21 : 30, : );
%         Z4 = s_E( 31 : 40, : );
%         Z5 = s_E( 41 : 50, : );
%         Z6 = s_E( 51 : 60, : );
%
%         M1 = mean( sum( Z1, 1 ), 2 );
%         M2 = mean( sum( Z2, 1 ), 2 );
%         M3 = mean( sum( Z3, 1 ), 2 );
%         M4 = mean( sum( Z4, 1 ), 2 );
%         M5 = mean( sum( Z5, 1 ), 2 );
%         M6 = mean( sum( Z6, 1 ), 2 );
%         [ ~, C1, H1 ] = patternMatching( S, Z1 );
%         [ ~, C2, H2 ] = patternMatching( S, Z2 );
%         [ ~, C3, H3 ] = patternMatching( S, Z3 );
%         [ ~, C4, H4 ] = patternMatching( S, Z4 );
%         [ ~, C5, H5 ] = patternMatching( S, Z5 );
%         [ ~, C6, H6 ] = patternMatching( S, Z6 );
%         C1 = C1( H1 > 0 );
%         C2 = C2( H2 > 0 );
%         C3 = C3( H3 > 0 );
%         C4 = C4( H4 > 0 );
%         C5 = C5( H5 > 0 );
%         C6 = C6( H6 > 0 );
%         if length( C1 ) >= 0.5 * length( H1 )
%             t1C = C1( 1 : end - 1 );
%             t2C = C1( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C1 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C1 = false;
%         end
%         if length( C2 ) >= 0.5 * length( H2 )
%             t1C = C2( 1 : end - 1 );
%             t2C = C2( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C2 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C2 = false;
%         end
%         if length( C3 ) >= 0.5 * length( H2 )
%             t1C = C3( 1 : end - 1 );
%             t2C = C3( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C3 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C3 = false;
%         end
%         if length( C4 ) >= 0.5 * length( H2 )
%             t1C = C4( 1 : end - 1 );
%             t2C = C4( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C4 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C4 = false;
%         end
%         if length( C5 ) >= 0.5 * length( H2 )
%             t1C = C5( 1 : end - 1 );
%             t2C = C5( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C5 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C5 = false;
%         end
%         if length( C6 ) >= 0.5 * length( H2 )
%             t1C = C6( 1 : end - 1 );
%             t2C = C6( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C6 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C6 = false;
%         end
%         H1 = mean( H1, 2 );
%         H2 = mean( H2, 2 );
%         H3 = mean( H3, 2 );
%         H4 = mean( H4, 2 );
%         H5 = mean( H5, 2 );
%         H6 = mean( H6, 2 );
%
%         M( ct_E, ct_I ) = mean( [ M1, M2, M3, M4, M5, M6 ], 2 );
%         C( ct_E, ct_I ) = all( [ C1, C2, C3, C4, C5, C6 ], 2 );
%         H( ct_E, ct_I ) = mean( [ H1, H2, H3, H4, H5, H6 ], 2 );
%
%         disp( [ 'E: ', num2str( ct_E ), ' / I: ', num2str( ct_I ) ] )
%     end; clear mI ct_I
% end; clear mE ct_E
%
% save( 'SNN_Simul6_1.mat', 'M', 'C', 'H', 's_E_all', 's_E_all_trans', 'pool_E', 'pool_I', '-v7.3' )
%

%% Result, SNN_Simul6
%
% load( 'SNN_Simul6_0.mat' )
% load( 'SNN_Simul6_1.mat' )
%
%
% figure( 'position', [ 100, 100, 300, 220 ] )
% imagesc( vars.J_EE, [ 0, 20 ] )
% colorbar% (pF)
% set( gca, 'xtick', [], 'ytick', [] )
%
% % ct_f = 0;
% % figure
% % for ct_E = 1 : 10
% %     for ct_I = 1 : 10
% %         ct_f = ct_f + 1;
% %         subplot( 10, 10, ct_f )
% %         imagesc( s_E_all_trans{ ct_E, ct_I } )
% %     end; clear ct_I
% % end; clear ct_E
%
%
% transMatU = zeros( 60, 3000 );
% ct_u = 0;
% for u = 1 : 60
%     transMatU( u, ct_u + [ 1 : 50 ] ) = 1;
%     ct_u = ct_u + 50;
% end; clear u ct_u
%
% spkT1 = [];% first spike timing
% for ct_E = 1 : length( pool_E )
%     for ct_I = 1 : length( pool_I )
%         s_E = s_E_all{ ct_E, ct_I };
%         s_E = ( transMatU * s_E );
%         spkT1( ct_E, ct_I ) = find( s_E( 1, : ) > 25, 1 );
%         clear s_E
%     end; clear ct_I
% end; clear ct_E
% spkT1 = mean( spkT1, 2 ) / 2;% (ms)
%
%
% colors = [ 1, 1, 1; 0.8, 0.8, 0.3 ];
%
% figure( 'position', [ 100, 100, 300, 200 ] )
% imagesc( 2 * ( 0.5 + 2 ) * pool_I, 1 * pool_E, C, [ -0.5, 1.5 ] )
% colormap( gca, colors )
% colorbar( 'ticks', [ 0, 1 ], 'ticklabels', { '', '' } )
% axis xy
% hold on
% ylabel( 'Excitatory rise param. (ms)' )
% xlabel( 'Approx. inhibition dur. (ms)' )
% title( 'States of parameters' )
%
% lw = 2;
% figure( 'position', [ 100, 100, 200, 200 ] )
% plot( pool_E, spkT1, '-', 'linewidth', lw )
% xlabel( 'Excitatory rise param. (ms)' )
% ylabel( 'Actual initiation dur. (ms)' )
%
% figure( 'position', [ 100, 100, 300, 200 ] )
% imagesc( 2 * ( 0.5 + 2 ) * pool_I, spkT1, C, [ -0.5, 1.5 ] )
% colormap( gca, colors )
% colorbar( 'ticks', [ 0, 1 ], 'ticklabels', { '', '' } )
% axis xy
% hold on
% ylabel( 'Actual initiation dur. (ms)' )
% xlabel( 'Approx. inhibition dur. (ms)' )
% title( 'States of parameters' )
%

%% SNN_Simul7_0
%
% clear vars
%
% vars.dt = 0.5;% time resolution (ms)
% vars.t_memory = 250;% time memory (ms)
%
% vars.tau_E = 20;% E neuron resting membrane time constant (ms)
% vars.tau_I = 20;% I neuron resting membrane time constant (ms)
% vars.E_L_E = -70;% E neuron resting potential (mV)
% vars.E_L_I = -62;% I neuron resting potential (mV)
% vars.Delta_T_E = 2;% E neuron EIF slope factor (mV)
% vars.Delta_T_I = 2;% I neuron EIF slope factor (mV)
% vars.C = 300;% Capacitance (pF)
% vars.E_E = 0;% E reverse potential (mV)
% vars.E_I = -75;% I reverse potential (mV)
% vars.V_T = -52;% threshold potential (mV)
% vars.A_T = 10;% Post spike threshold potential increase (mV)
% vars.tau_T = 30;% Adaptive threshold time scale (ms)
% vars.V_re = -60;% Reset potential (mV)
% vars.tau_abs = 1;% Absolute refractory period (ms)
% vars.a_w = 4;% Subthreshold adaptation (nS)
% vars.b_w = 0.805;% Spike-triggered adaptation (pA)
% vars.tau_w = 150;% Spike-triggered adaptation time scale (ms)
% vars.V_ap = 20;% action potential threshold potention (mV)
% vars.V_lb = -75;% lower bound potential (mV) := vars.E_I
%
% vars.N_E = [ 500, 500, 500, 500, 500, 500, 500 ];% the number of excitatory neurons
% vars.N_I = [ 125, 125, 125, 125, 125, 125, 125 ];% the number of excitatory neurons
% % vars.N_E = [ 1000, 1000 ];% the number of excitatory neurons
% % vars.N_I = [ 250, 250 ];% the number of excitatory neurons
% % vars.p_C = 2.5 * 0.2;% connection probability via chemical synapses
% vars.p_C = 0.2;% connection probability via chemical synapses
% vars.scaleWeight = 0.04;% parameter for re-scaling of synaptic weights
% % vars.scaleWeight = 0.12;% parameter for re-scaling of synaptic weights
% % vars.scaleWeight = 0.15;% parameter for re-scaling of synaptic weights
% % vars.scaleWeight = 0.2;% parameter for re-scaling of synaptic weights
% % vars.p_E = 0.2;% connection probability via electrical synapses
% % vars.tau_r_E = 1;% rise time for E synapses (ms)
% % vars.tau_d_E = 6;% decay time for E synapses (ms)
% % vars.tau_r_I = 0.5;% rise time for I synapses (ms)
% % vars.tau_d_I = 2;% decay time for I synapses (ms)
% vars.r_ext_EE = 12 * 500;% Rate of external input to E neurons (Hz)
% vars.r_ext_IE = 0 * 500;% Rate of external input to I neurons (Hz)
% vars.J_EE_min = 1.78;% minimum E to E synaptic weight (pF)
% vars.J_EE_max = 21.4;% maximum E to E synaptic weight (pF)
% vars.J_EE_O = 2.76;% initial E to E synaptic weight (pF)
% vars.J_EE_O1 = 2.76 * 8;% initial E to E synaptic weight (pF)
% vars.J_EI_min = 48.7;% minimum I to E synaptic weight (pF)
% vars.J_EI_max = 243;% maximum I to E synaptic weight (pF)
% vars.J_EI_O = 3 * 48.7;% initial I to E synaptic weight (pF)
% vars.J_IE_val = 1.27;% value of E to I synaptic weight (pF)
% vars.J_II_val = 16.2;% value of I to I synaptic weight (pF)
% % vars.J_GII_val = 0;% value of I to I gap junction weight (pF)
% vars.J_ext_EE = vars.J_EE_min;
% vars.J_ext_IE = vars.J_IE_val;
% vars.period_scaling = 20;% period of synaptic scaling (ms)
%
% vars.A_LTD = 0.0008;% Long-term depression (LTD) strength (pA mV^-1)
% vars.A_LTP = 0.0014;% Long-term potentiation (LTP) strength (pA mV^-2)
% vars.theta_LTD = -70;% Threshold to recruit LTD (mV)
% vars.theta_LTP = -49;% Threshold to recruit LTP (mV)
% vars.tau_u = 10;% Time constant of low-pass filtered membrane voltage (for LTD) (ms)
% vars.tau_v = 7;% Time constant of low-pass filtered membrane voltage (for LTP) (ms)
% vars.tau_x = 15;% Time constant low-pass filtered spike train (for LTP) (ms)
%
% vars.tau_y = 20;% Time constant of low-pass filtered spike train (ms)
% vars.eta = 1;% Synaptic plasticity learning rate (pA)
% vars.r_O = 3;% Target firing rate (Hz)
%
% % vars.triplet_tau_plus = 16.8;% parameter of triplet plasticity (ms)
% % vars.triplet_tau_minus = 33.7;% parameter of triplet plasticity (ms)
% % vars.triplet_tau_x = 101;% parameter of triplet plasticity (ms)
% % vars.triplet_tau_y = 125;% parameter of triplet plasticity (ms)
% % vars.triplet_A_2_plus = 7.5 * 10 ^ (-10);% parameter of triplet plasticity
% % vars.triplet_A_3_plus = 9.3 * 10 ^ (-3);% parameter of triplet plasticity
% % vars.triplet_A_2_minus = 7 * 10 ^ (-3);% parameter of triplet plasticity
% % vars.triplet_A_3_minus = 2.3 * 10 ^ (-4);% parameter of triplet plasticity
%
% vars.V_E = [];% voltages of E (mV)
% vars.V_I = [];% voltages of I (mV)
% vars.Vu_E = [];% low-pass filtered valtage of E (mV)
% vars.Vv_E = [];% low-pass filtered valtage of E (mV)
%
% vars.s_E = [];% spike train of E
% vars.sx_E = [];% low-pass filtered spike train of E
% vars.sy_E = [];% low-pass filtered spike train of E
% vars.s_I = [];% spike train of I
% vars.sy_I = [];% low-pass filtered spike train of I
%
% vars.V_T_E = [];% thereshold of E
% vars.ww_E = [];% adaptation current of E
%
% vars.s_ext_EE = [];% external stimulation to E
% vars.s_ext_IE = [];% external stimulation to I
%
% vars.ct_abs_E = [];% count of absolute refractory period for E neurons
% vars.ct_abs_I = [];% count of absolute refractory period for I neurons
%
% % vars.triplet_r_1 = [];% variable of triplet plasticity
% % vars.triplet_r_2 = [];% variable of triplet plasticity
% % vars.triplet_o_1 = [];% variable of triplet plasticity
% % vars.triplet_o_2 = [];% variable of triplet plasticity
%
% vars.connScale = [ 1, 1, 7 ];% [ orthodromic, antidromic, within area ]
% vars.connects = zeros( length( vars.N_E ), length( vars.N_E ) );% vars.connects( j, i ): from area i to area j.
% % vars.connects( 2, 1 ) = vars.connScale( 1 ) * vars.N_E( 1 );
% % vars.connects( 1, 2 ) = vars.connScale( 2 ) * vars.N_E( 2 );
% % vars.connects( 1, 1 ) = vars.connScale( 3 ) * vars.N_E( 1 );
% % vars.connects( 2, 2 ) = vars.connScale( 3 ) * vars.N_E( 2 );
% % vars.connects = vars.connects ./ sum( vars.connects, 2 );
% % vars.connects( isnan( vars.connects ) ) = 0;
% % vars.connects( vars.connects > 0 ) = 1;
% vars.connects( 2, 1 ) = vars.connScale( 1 );
% vars.connects( 3, 2 ) = vars.connScale( 1 );
% vars.connects( 4, 3 ) = vars.connScale( 1 );
% vars.connects( 5, 4 ) = vars.connScale( 1 );
% vars.connects( 6, 5 ) = vars.connScale( 1 );
% vars.connects( 7, 6 ) = vars.connScale( 1 );
% vars.connects( 1, 7 ) = vars.connScale( 2 );
% vars.connects( 1, 1 ) = vars.connScale( 3 );
% vars.connects( 2, 2 ) = vars.connScale( 3 );
% vars.connects( 3, 3 ) = vars.connScale( 3 );
% vars.connects( 4, 4 ) = vars.connScale( 3 );
% vars.connects( 5, 5 ) = vars.connScale( 3 );
% vars.connects( 6, 6 ) = vars.connScale( 3 );
% vars.connects( 7, 7 ) = vars.connScale( 3 );
%
% vars.tau_r_E = 1 * 2;% rise time for E synapses (ms)
% vars.tau_d_E = 6 * 2;% decay time for E synapses (ms)
% vars.tau_r_I = 0.5 * 10;% rise time for I synapses (ms)
% vars.tau_d_I = 2 * 10;% decay time for I synapses (ms)
% vars.nChannels = 10;
% vars = SNN_initiation0( vars );
%
%
% save( 'SNN_Simul7_0.mat', 'vars' )
% clear vars
%
%
% % dt = vars.dt;
% % t_memory = vars.t_memory;
% % t_memory_L = t_memory / dt;
% % tau_r_E = vars.tau_r_E;
% % tau_d_E = vars.tau_d_E;
% % tau_r_I = vars.tau_r_I;
% % tau_d_I = vars.tau_d_I;
% % t_domain_memory = 0 : dt : t_memory - dt;
% % F_E = ( 1 / ( tau_d_E - tau_r_E ) ) * ( exp( -t_domain_memory / tau_d_E ) - exp( -t_domain_memory / tau_r_E ) );
% % F_I = ( 1 / ( tau_d_I - tau_r_I ) ) * ( exp( -t_domain_memory / tau_d_I ) - exp( -t_domain_memory / tau_r_I ) );
% % figure; hold on; plot(F_E); plot(F_I)
%

%% SNN_Simul7_1
%
% transMatU = zeros( 70, 3500 );
% ct_u = 0;
% for u = 1 : 70
%     transMatU( u, ct_u + [ 1 : 50 ] ) = 1;
%     ct_u = ct_u + 50;
% end; clear u ct_u
%
% transMatT = zeros( 9600, 480 );
% ct_t = 0;
% for t = 1 : 480
%     transMatT( ct_t + [ 1 : 20 ], t ) = 1;
%     ct_t = ct_t + 20;
% end; clear t ct_t
%
% n_stim = 10;
% S = eye( n_stim );
%
%
% pool_E = [ 7 : 2 : 25 ];
% pool_I = [ 3 : 3 : 30 ];
%
% M = [];
% C = [];
% H = [];
% s_E_all = {};
% s_E_all_trans = {};
%
% ct_E = 0;
% for mE = pool_E
%     ct_E = ct_E + 1;
%     ct_I = 0;
%     for mI = pool_I
%         ct_I = ct_I + 1;
%
%         load( 'SNN_Simul7_0.mat' )
%
%         vars.tau_r_E = 1 * mE;% rise time for E synapses (ms)
%         vars.tau_d_E = 6 * mE * 0.3;% decay time for E synapses (ms)
%         vars.tau_r_I = 0.5 * mI;% rise time for I synapses (ms)
%         vars.tau_d_I = 2 * mI;% decay time for I synapses (ms)
%
%         vars.mode = 0.5;% the current balance between synaptic weights ( 1:= lateral, 0:= recurrent ), \in [ 0, 1 ].
%         vars.infoStims = [ 4800, 120, 40, 10, 1 ];% recording duration (ms), stimulation duration (ms), stimulation delay (ms), the number of assemblies, the number of stimulations ];
%         vars.inputs = [ 1, 2, 3, 4, 5, 6, 7 ];% in area
%         vars = SNN_stimulation( vars, 0 );
%         [ vars, records ] = SNN_simulation_GPU( vars, 'plasticity_off', 'records_on' );
%
%         s_E = records.s_E;
%         s_E_all{ ct_E, ct_I } = s_E;
%         s_E = ( transMatU * s_E ) * transMatT;
%         s_E_all_trans{ ct_E, ct_I } = s_E;
%
%         % figure
%         % imagesc( vars.ext_stim_E ); colorbar
%         % figure
%         % imagesc( records.s_E ); colorbar
%         % figure
%         % imagesc( s_E ); colorbar
%
%
%         Z1 = s_E( 1 : 10, : );
%         Z2 = s_E( 11 : 20, : );
%         Z3 = s_E( 21 : 30, : );
%         Z4 = s_E( 31 : 40, : );
%         Z5 = s_E( 41 : 50, : );
%         Z6 = s_E( 51 : 60, : );
%         Z7 = s_E( 61 : 70, : );
%
%         M1 = mean( sum( Z1, 1 ), 2 );
%         M2 = mean( sum( Z2, 1 ), 2 );
%         M3 = mean( sum( Z3, 1 ), 2 );
%         M4 = mean( sum( Z4, 1 ), 2 );
%         M5 = mean( sum( Z5, 1 ), 2 );
%         M6 = mean( sum( Z6, 1 ), 2 );
%         M7 = mean( sum( Z7, 1 ), 2 );
%         [ ~, C1, H1 ] = patternMatching( S, Z1 );
%         [ ~, C2, H2 ] = patternMatching( S, Z2 );
%         [ ~, C3, H3 ] = patternMatching( S, Z3 );
%         [ ~, C4, H4 ] = patternMatching( S, Z4 );
%         [ ~, C5, H5 ] = patternMatching( S, Z5 );
%         [ ~, C6, H6 ] = patternMatching( S, Z6 );
%         [ ~, C7, H7 ] = patternMatching( S, Z7 );
%         C1 = C1( H1 > 0 );
%         C2 = C2( H2 > 0 );
%         C3 = C3( H3 > 0 );
%         C4 = C4( H4 > 0 );
%         C5 = C5( H5 > 0 );
%         C6 = C6( H6 > 0 );
%         C7 = C7( H7 > 0 );
%         if length( C1 ) >= 0.5 * length( H1 )
%             t1C = C1( 1 : end - 1 );
%             t2C = C1( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C1 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C1 = false;
%         end
%         if length( C2 ) >= 0.5 * length( H2 )
%             t1C = C2( 1 : end - 1 );
%             t2C = C2( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C2 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C2 = false;
%         end
%         if length( C3 ) >= 0.5 * length( H2 )
%             t1C = C3( 1 : end - 1 );
%             t2C = C3( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C3 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C3 = false;
%         end
%         if length( C4 ) >= 0.5 * length( H2 )
%             t1C = C4( 1 : end - 1 );
%             t2C = C4( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C4 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C4 = false;
%         end
%         if length( C5 ) >= 0.5 * length( H2 )
%             t1C = C5( 1 : end - 1 );
%             t2C = C5( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C5 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C5 = false;
%         end
%         if length( C6 ) >= 0.5 * length( H2 )
%             t1C = C6( 1 : end - 1 );
%             t2C = C6( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C6 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C6 = false;
%         end
%         if length( C7 ) >= 0.5 * length( H2 )
%             t1C = C7( 1 : end - 1 );
%             t2C = C7( 2 : end );
%             idxC = [];
%             for n = 1 : n_stim
%                 if n < n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%                 elseif n == n_stim
%                     idx = t1C == n;
%                     idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%                 end
%                 if ~isempty( idx )
%                     idxC = [ idxC, all( idx, 2 ) ];
%                 else
%                     idxC = [ idxC, false ];
%                 end
%                 clear idx
%             end; clear n
%             C7 = mean( idxC, 2 ) >= 0.5;
%             clear t1C t2C idxC
%         else
%             C7 = false;
%         end
%         H1 = mean( H1, 2 );
%         H2 = mean( H2, 2 );
%         H3 = mean( H3, 2 );
%         H4 = mean( H4, 2 );
%         H5 = mean( H5, 2 );
%         H6 = mean( H6, 2 );
%         H7 = mean( H7, 2 );
%
%         M( ct_E, ct_I ) = mean( [ M1, M2, M3, M4, M5, M6, M7 ], 2 );
%         C( ct_E, ct_I ) = all( [ C1, C2, C3, C4, C5, C6, C7 ], 2 );
%         H( ct_E, ct_I ) = mean( [ H1, H2, H3, H4, H5, H6, H7 ], 2 );
%
%         disp( [ 'E: ', num2str( ct_E ), ' / I: ', num2str( ct_I ) ] )
%     end; clear mI ct_I
% end; clear mE ct_E
%
% save( 'SNN_Simul7_1.mat', 'M', 'C', 'H', 's_E_all', 's_E_all_trans', 'pool_E', 'pool_I', '-v7.3' )
%

%% Result, SNN_Simul7
%
% load( 'SNN_Simul7_0.mat' )
% load( 'SNN_Simul7_1.mat' )
%
%
% figure( 'position', [ 100, 100, 300, 220 ] )
% imagesc( vars.J_EE, [ 0, 20 ] )
% colorbar% (pF)
% set( gca, 'xtick', [], 'ytick', [] )
%
% % ct_f = 0;
% % figure
% % for ct_E = 1 : 10
% %     for ct_I = 1 : 10
% %         ct_f = ct_f + 1;
% %         subplot( 10, 10, ct_f )
% %         imagesc( s_E_all_trans{ ct_E, ct_I } )
% %     end; clear ct_I
% % end; clear ct_E
%
%
% transMatU = zeros( 70, 3500 );
% ct_u = 0;
% for u = 1 : 70
%     transMatU( u, ct_u + [ 1 : 50 ] ) = 1;
%     ct_u = ct_u + 50;
% end; clear u ct_u
%
% spkT1 = [];% first spike timing
% for ct_E = 1 : length( pool_E )
%     for ct_I = 1 : length( pool_I )
%         s_E = s_E_all{ ct_E, ct_I };
%         s_E = ( transMatU * s_E );
%         spkT1( ct_E, ct_I ) = find( s_E( 1, : ) > 25, 1 );
%         clear s_E
%     end; clear ct_I
% end; clear ct_E
% spkT1 = mean( spkT1, 2 ) / 2;% (ms)
%
%
% colors = [ 1, 1, 1; 0.8, 0.8, 0.3 ];
%
% figure( 'position', [ 100, 100, 300, 200 ] )
% imagesc( 2 * ( 0.5 + 2 ) * pool_I, 1 * pool_E, C, [ -0.5, 1.5 ] )
% colormap( gca, colors )
% colorbar( 'ticks', [ 0, 1 ], 'ticklabels', { '', '' } )
% axis xy
% hold on
% ylabel( 'Excitatory rise param. (ms)' )
% xlabel( 'Approx. inhibition dur. (ms)' )
% title( 'States of parameters' )
%
% lw = 2;
% figure( 'position', [ 100, 100, 200, 200 ] )
% plot( pool_E, spkT1, '-', 'linewidth', lw )
% xlabel( 'Excitatory rise param. (ms)' )
% ylabel( 'Actual initiation dur. (ms)' )
%
% figure( 'position', [ 100, 100, 300, 200 ] )
% imagesc( 2 * ( 0.5 + 2 ) * pool_I, spkT1, C, [ -0.5, 1.5 ] )
% colormap( gca, colors )
% colorbar( 'ticks', [ 0, 1 ], 'ticklabels', { '', '' } )
% axis xy
% hold on
% ylabel( 'Actual initiation dur. (ms)' )
% xlabel( 'Approx. inhibition dur. (ms)' )
% title( 'States of parameters' )
%

%% SNN_Simul2_2
%
% transMatU = zeros( 20, 1000 );
% ct_u = 0;
% for u = 1 : 20
%     transMatU( u, ct_u + [ 1 : 50 ] ) = 1;
%     ct_u = ct_u + 50;
% end; clear u ct_u
%
% transMatT = zeros( 9600, 480 );
% ct_t = 0;
% for t = 1 : 480
%     transMatT( ct_t + [ 1 : 20 ], t ) = 1;
%     ct_t = ct_t + 20;
% end; clear t ct_t
%
% n_stim = 10;
% S = eye( n_stim );
%
%
% pool_E = [ 7 : 2 : 25 ];
% pool_I = [ 3 : 3 : 30 ];
%
% M = [];
% C = [];
% H = [];
% s_E_all = {};
% s_E_all_trans = {};
%
% mE = pool_E( 5 );
% mI = pool_I( 7 );
% ct_noiseLevel = 0;
% for noiseLevel = [ 0 : 0.2 : 4 ]
%
%     ct_noiseLevel = ct_noiseLevel + 1;
%
%     load( 'SNN_Simul2_0.mat' )
%
%     vars.tau_r_E = 1 * mE;% rise time for E synapses (ms)
%     vars.tau_d_E = 6 * mE * 0.3;% decay time for E synapses (ms)
%     vars.tau_r_I = 0.5 * mI;% rise time for I synapses (ms)
%     vars.tau_d_I = 2 * mI;% decay time for I synapses (ms)
%
%     vars.mode = 0.5;% the current balance between synaptic weights ( 1:= lateral, 0:= recurrent ), \in [ 0, 1 ].
%     vars.infoStims = [ 4800, 120, 40, 10, 1 ];% recording duration (ms), stimulation duration (ms), stimulation delay (ms), the number of assemblies, the number of stimulations ];
%     vars.inputs = [ 1, 2 ];% in area
%     vars = SNN_stimulation( vars, noiseLevel );
%     [ vars, records ] = SNN_simulation( vars, 'plasticity_off', 'records_on' );
%
%     s_E = records.s_E;
%     s_E_all{ ct_noiseLevel, 1 } = s_E;
%     s_E = ( transMatU * s_E ) * transMatT;
%     s_E_all_trans{ ct_noiseLevel, 1 } = s_E;
%
%     % figure
%     % imagesc( vars.ext_stim_E ); colorbar
%     % figure
%     % imagesc( records.s_E ); colorbar
%     % figure
%     % imagesc( s_E ); colorbar
%
%
%     Z1 = s_E( 1 : 10, : );
%     Z2 = s_E( 11 : 20, : );
%
%     M1 = mean( sum( Z1, 1 ), 2 );
%     M2 = mean( sum( Z2, 1 ), 2 );
%     [ ~, C1, H1 ] = patternMatching( S, Z1 );
%     [ ~, C2, H2 ] = patternMatching( S, Z2 );
%     C1 = C1( H1 > 0 );
%     C2 = C2( H2 > 0 );
%     if length( C1 ) >= 0.5 * length( H1 )
%         t1C = C1( 1 : end - 1 );
%         t2C = C1( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C1 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C1 = false;
%     end
%     if length( C2 ) >= 0.5 * length( H2 )
%         t1C = C2( 1 : end - 1 );
%         t2C = C2( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C2 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C2 = false;
%     end
%     H1 = mean( H1, 2 );
%     H2 = mean( H2, 2 );
%
%     M( ct_noiseLevel, 1 ) = mean( [ M1, M2 ], 2 );
%     C( ct_noiseLevel, 1 ) = all( [ C1, C2 ], 2 );
%     H( ct_noiseLevel, 1 ) = mean( [ H1, H2 ], 2 );
%
%     disp( [ 'NoiseLevel: ', num2str( noiseLevel ) ] )
% end; clear noiseLevel
%
% save( 'SNN_Simul2_2.mat', 'M', 'C', 'H', 's_E_all', 's_E_all_trans', 'pool_E', 'pool_I', '-v7.3' )
%

%% SNN_Simul3_2
%
% transMatU = zeros( 30, 1500 );
% ct_u = 0;
% for u = 1 : 30
%     transMatU( u, ct_u + [ 1 : 50 ] ) = 1;
%     ct_u = ct_u + 50;
% end; clear u ct_u
%
% transMatT = zeros( 9600, 480 );
% ct_t = 0;
% for t = 1 : 480
%     transMatT( ct_t + [ 1 : 20 ], t ) = 1;
%     ct_t = ct_t + 20;
% end; clear t ct_t
%
% n_stim = 10;
% S = eye( n_stim );
%
%
% pool_E = [ 7 : 2 : 25 ];
% pool_I = [ 3 : 3 : 30 ];
%
% M = [];
% C = [];
% H = [];
% s_E_all = {};
% s_E_all_trans = {};
%
% mE = pool_E( 5 );
% mI = pool_I( 7 );
% ct_noiseLevel = 0;
% for noiseLevel = [ 0 : 0.2 : 4 ]
%
%     ct_noiseLevel = ct_noiseLevel + 1;
%
%     load( 'SNN_Simul3_0.mat' )
%
%     vars.tau_r_E = 1 * mE;% rise time for E synapses (ms)
%     vars.tau_d_E = 6 * mE * 0.3;% decay time for E synapses (ms)
%     vars.tau_r_I = 0.5 * mI;% rise time for I synapses (ms)
%     vars.tau_d_I = 2 * mI;% decay time for I synapses (ms)
%
%     vars.mode = 0.5;% the current balance between synaptic weights ( 1:= lateral, 0:= recurrent ), \in [ 0, 1 ].
%     vars.infoStims = [ 4800, 120, 40, 10, 1 ];% recording duration (ms), stimulation duration (ms), stimulation delay (ms), the number of assemblies, the number of stimulations ];
%     vars.inputs = [ 1, 2, 3 ];% in area
%     vars = SNN_stimulation( vars, noiseLevel );
%     [ vars, records ] = SNN_simulation( vars, 'plasticity_off', 'records_on' );
%
%     s_E = records.s_E;
%     s_E_all{ ct_noiseLevel, 1 } = s_E;
%     s_E = ( transMatU * s_E ) * transMatT;
%     s_E_all_trans{ ct_noiseLevel, 1 } = s_E;
%
%     % figure
%     % imagesc( vars.ext_stim_E ); colorbar
%     % figure
%     % imagesc( records.s_E ); colorbar
%     % figure
%     % imagesc( s_E ); colorbar
%
%
%     Z1 = s_E( 1 : 10, : );
%     Z2 = s_E( 11 : 20, : );
%     Z3 = s_E( 21 : 30, : );
%
%     M1 = mean( sum( Z1, 1 ), 2 );
%     M2 = mean( sum( Z2, 1 ), 2 );
%     M3 = mean( sum( Z3, 1 ), 2 );
%     [ ~, C1, H1 ] = patternMatching( S, Z1 );
%     [ ~, C2, H2 ] = patternMatching( S, Z2 );
%     [ ~, C3, H3 ] = patternMatching( S, Z3 );
%     C1 = C1( H1 > 0 );
%     C2 = C2( H2 > 0 );
%     C3 = C3( H3 > 0 );
%     if length( C1 ) >= 0.5 * length( H1 )
%         t1C = C1( 1 : end - 1 );
%         t2C = C1( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C1 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C1 = false;
%     end
%     if length( C2 ) >= 0.5 * length( H2 )
%         t1C = C2( 1 : end - 1 );
%         t2C = C2( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C2 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C2 = false;
%     end
%     if length( C3 ) >= 0.5 * length( H2 )
%         t1C = C3( 1 : end - 1 );
%         t2C = C3( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C3 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C3 = false;
%     end
%     H1 = mean( H1, 2 );
%     H2 = mean( H2, 2 );
%     H3 = mean( H3, 2 );
%
%     M( ct_noiseLevel, 1 ) = mean( [ M1, M2, M3 ], 2 );
%     C( ct_noiseLevel, 1 ) = all( [ C1, C2, C3 ], 2 );
%     H( ct_noiseLevel, 1 ) = mean( [ H1, H2, H3 ], 2 );
%
%     disp( [ 'NoiseLevel: ', num2str( noiseLevel ) ] )
% end; clear noiseLevel
%
% save( 'SNN_Simul3_2.mat', 'M', 'C', 'H', 's_E_all', 's_E_all_trans', 'pool_E', 'pool_I', '-v7.3' )
%

%% SNN_Simul4_2
%
% transMatU = zeros( 40, 2000 );
% ct_u = 0;
% for u = 1 : 40
%     transMatU( u, ct_u + [ 1 : 50 ] ) = 1;
%     ct_u = ct_u + 50;
% end; clear u ct_u
%
% transMatT = zeros( 9600, 480 );
% ct_t = 0;
% for t = 1 : 480
%     transMatT( ct_t + [ 1 : 20 ], t ) = 1;
%     ct_t = ct_t + 20;
% end; clear t ct_t
%
% n_stim = 10;
% S = eye( n_stim );
%
%
% pool_E = [ 7 : 2 : 25 ];
% pool_I = [ 3 : 3 : 30 ];
%
% M = [];
% C = [];
% H = [];
% s_E_all = {};
% s_E_all_trans = {};
%
% mE = pool_E( 5 );
% mI = pool_I( 7 );
% ct_noiseLevel = 0;
% for noiseLevel = [ 0 : 0.2 : 4 ]
%
%     ct_noiseLevel = ct_noiseLevel + 1;
%
%     load( 'SNN_Simul4_0.mat' )
%
%     vars.tau_r_E = 1 * mE;% rise time for E synapses (ms)
%     vars.tau_d_E = 6 * mE * 0.3;% decay time for E synapses (ms)
%     vars.tau_r_I = 0.5 * mI;% rise time for I synapses (ms)
%     vars.tau_d_I = 2 * mI;% decay time for I synapses (ms)
%
%     vars.mode = 0.5;% the current balance between synaptic weights ( 1:= lateral, 0:= recurrent ), \in [ 0, 1 ].
%     vars.infoStims = [ 4800, 120, 40, 10, 1 ];% recording duration (ms), stimulation duration (ms), stimulation delay (ms), the number of assemblies, the number of stimulations ];
%     vars.inputs = [ 1, 2, 3, 4 ];% in area
%     vars = SNN_stimulation( vars, noiseLevel );
%     [ vars, records ] = SNN_simulation( vars, 'plasticity_off', 'records_on' );
%
%     s_E = records.s_E;
%     s_E_all{ ct_noiseLevel, 1 } = s_E;
%     s_E = ( transMatU * s_E ) * transMatT;
%     s_E_all_trans{ ct_noiseLevel, 1 } = s_E;
%
%     % figure
%     % imagesc( vars.ext_stim_E ); colorbar
%     % figure
%     % imagesc( records.s_E ); colorbar
%     % figure
%     % imagesc( s_E ); colorbar
%
%
%     Z1 = s_E( 1 : 10, : );
%     Z2 = s_E( 11 : 20, : );
%     Z3 = s_E( 21 : 30, : );
%     Z4 = s_E( 31 : 40, : );
%
%     M1 = mean( sum( Z1, 1 ), 2 );
%     M2 = mean( sum( Z2, 1 ), 2 );
%     M3 = mean( sum( Z3, 1 ), 2 );
%     M4 = mean( sum( Z4, 1 ), 2 );
%     [ ~, C1, H1 ] = patternMatching( S, Z1 );
%     [ ~, C2, H2 ] = patternMatching( S, Z2 );
%     [ ~, C3, H3 ] = patternMatching( S, Z3 );
%     [ ~, C4, H4 ] = patternMatching( S, Z4 );
%     C1 = C1( H1 > 0 );
%     C2 = C2( H2 > 0 );
%     C3 = C3( H3 > 0 );
%     C4 = C4( H4 > 0 );
%     if length( C1 ) >= 0.5 * length( H1 )
%         t1C = C1( 1 : end - 1 );
%         t2C = C1( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C1 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C1 = false;
%     end
%     if length( C2 ) >= 0.5 * length( H2 )
%         t1C = C2( 1 : end - 1 );
%         t2C = C2( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C2 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C2 = false;
%     end
%     if length( C3 ) >= 0.5 * length( H2 )
%         t1C = C3( 1 : end - 1 );
%         t2C = C3( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C3 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C3 = false;
%     end
%     if length( C4 ) >= 0.5 * length( H2 )
%         t1C = C4( 1 : end - 1 );
%         t2C = C4( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C4 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C4 = false;
%     end
%     H1 = mean( H1, 2 );
%     H2 = mean( H2, 2 );
%     H3 = mean( H3, 2 );
%     H4 = mean( H4, 2 );
%
%     M( ct_noiseLevel, 1 ) = mean( [ M1, M2, M3, M4 ], 2 );
%     C( ct_noiseLevel, 1 ) = all( [ C1, C2, C3, C4 ], 2 );
%     H( ct_noiseLevel, 1 ) = mean( [ H1, H2, H3, H4 ], 2 );
%
%     disp( [ 'NoiseLevel: ', num2str( noiseLevel ) ] )
% end; clear noiseLevel
%
% save( 'SNN_Simul4_2.mat', 'M', 'C', 'H', 's_E_all', 's_E_all_trans', 'pool_E', 'pool_I', '-v7.3' )
%

%% SNN_Simul5_2
%
% transMatU = zeros( 50, 2500 );
% ct_u = 0;
% for u = 1 : 50
%     transMatU( u, ct_u + [ 1 : 50 ] ) = 1;
%     ct_u = ct_u + 50;
% end; clear u ct_u
%
% transMatT = zeros( 9600, 480 );
% ct_t = 0;
% for t = 1 : 480
%     transMatT( ct_t + [ 1 : 20 ], t ) = 1;
%     ct_t = ct_t + 20;
% end; clear t ct_t
%
% n_stim = 10;
% S = eye( n_stim );
%
%
% pool_E = [ 7 : 2 : 25 ];
% pool_I = [ 3 : 3 : 30 ];
%
% M = [];
% C = [];
% H = [];
% s_E_all = {};
% s_E_all_trans = {};
%
% mE = pool_E( 5 );
% mI = pool_I( 7 );
% ct_noiseLevel = 0;
% for noiseLevel = [ 0 : 0.2 : 4 ]
%
%     ct_noiseLevel = ct_noiseLevel + 1;
%
%     load( 'SNN_Simul5_0.mat' )
%
%     vars.tau_r_E = 1 * mE;% rise time for E synapses (ms)
%     vars.tau_d_E = 6 * mE * 0.3;% decay time for E synapses (ms)
%     vars.tau_r_I = 0.5 * mI;% rise time for I synapses (ms)
%     vars.tau_d_I = 2 * mI;% decay time for I synapses (ms)
%
%     vars.mode = 0.5;% the current balance between synaptic weights ( 1:= lateral, 0:= recurrent ), \in [ 0, 1 ].
%     vars.infoStims = [ 4800, 120, 40, 10, 1 ];% recording duration (ms), stimulation duration (ms), stimulation delay (ms), the number of assemblies, the number of stimulations ];
%     vars.inputs = [ 1, 2, 3, 4, 5 ];% in area
%     vars = SNN_stimulation( vars, noiseLevel );
%     [ vars, records ] = SNN_simulation_GPU( vars, 'plasticity_off', 'records_on' );
%
%     s_E = records.s_E;
%     s_E_all{ ct_noiseLevel, 1 } = s_E;
%     s_E = ( transMatU * s_E ) * transMatT;
%     s_E_all_trans{ ct_noiseLevel, 1 } = s_E;
%
%     % figure
%     % imagesc( vars.ext_stim_E ); colorbar
%     % figure
%     % imagesc( records.s_E ); colorbar
%     % figure
%     % imagesc( s_E ); colorbar
%
%
%     Z1 = s_E( 1 : 10, : );
%     Z2 = s_E( 11 : 20, : );
%     Z3 = s_E( 21 : 30, : );
%     Z4 = s_E( 31 : 40, : );
%     Z5 = s_E( 41 : 50, : );
%
%     M1 = mean( sum( Z1, 1 ), 2 );
%     M2 = mean( sum( Z2, 1 ), 2 );
%     M3 = mean( sum( Z3, 1 ), 2 );
%     M4 = mean( sum( Z4, 1 ), 2 );
%     M5 = mean( sum( Z5, 1 ), 2 );
%     [ ~, C1, H1 ] = patternMatching( S, Z1 );
%     [ ~, C2, H2 ] = patternMatching( S, Z2 );
%     [ ~, C3, H3 ] = patternMatching( S, Z3 );
%     [ ~, C4, H4 ] = patternMatching( S, Z4 );
%     [ ~, C5, H5 ] = patternMatching( S, Z5 );
%     C1 = C1( H1 > 0 );
%     C2 = C2( H2 > 0 );
%     C3 = C3( H3 > 0 );
%     C4 = C4( H4 > 0 );
%     C5 = C5( H5 > 0 );
%     if length( C1 ) >= 0.5 * length( H1 )
%         t1C = C1( 1 : end - 1 );
%         t2C = C1( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C1 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C1 = false;
%     end
%     if length( C2 ) >= 0.5 * length( H2 )
%         t1C = C2( 1 : end - 1 );
%         t2C = C2( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C2 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C2 = false;
%     end
%     if length( C3 ) >= 0.5 * length( H2 )
%         t1C = C3( 1 : end - 1 );
%         t2C = C3( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C3 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C3 = false;
%     end
%     if length( C4 ) >= 0.5 * length( H2 )
%         t1C = C4( 1 : end - 1 );
%         t2C = C4( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C4 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C4 = false;
%     end
%     if length( C5 ) >= 0.5 * length( H2 )
%         t1C = C5( 1 : end - 1 );
%         t2C = C5( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C5 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C5 = false;
%     end
%     H1 = mean( H1, 2 );
%     H2 = mean( H2, 2 );
%     H3 = mean( H3, 2 );
%     H4 = mean( H4, 2 );
%     H5 = mean( H5, 2 );
%
%     M( ct_noiseLevel, 1 ) = mean( [ M1, M2, M3, M4, M5 ], 2 );
%     C( ct_noiseLevel, 1 ) = all( [ C1, C2, C3, C4, C5 ], 2 );
%     H( ct_noiseLevel, 1 ) = mean( [ H1, H2, H3, H4, H5 ], 2 );
%
%     disp( [ 'NoiseLevel: ', num2str( noiseLevel ) ] )
% end; clear noiseLevel
%
% save( 'SNN_Simul5_2.mat', 'M', 'C', 'H', 's_E_all', 's_E_all_trans', 'pool_E', 'pool_I', '-v7.3' )
%

%% SNN_Simul6_2
%
% transMatU = zeros( 60, 3000 );
% ct_u = 0;
% for u = 1 : 60
%     transMatU( u, ct_u + [ 1 : 50 ] ) = 1;
%     ct_u = ct_u + 50;
% end; clear u ct_u
%
% transMatT = zeros( 9600, 480 );
% ct_t = 0;
% for t = 1 : 480
%     transMatT( ct_t + [ 1 : 20 ], t ) = 1;
%     ct_t = ct_t + 20;
% end; clear t ct_t
%
% n_stim = 10;
% S = eye( n_stim );
%
%
% pool_E = [ 7 : 2 : 25 ];
% pool_I = [ 3 : 3 : 30 ];
%
% M = [];
% C = [];
% H = [];
% s_E_all = {};
% s_E_all_trans = {};
%
% mE = pool_E( 5 );
% mI = pool_I( 7 );
% ct_noiseLevel = 0;
% for noiseLevel = [ 0 : 0.2 : 4 ]
%
%     ct_noiseLevel = ct_noiseLevel + 1;
%
%     load( 'SNN_Simul6_0.mat' )
%
%     vars.tau_r_E = 1 * mE;% rise time for E synapses (ms)
%     vars.tau_d_E = 6 * mE * 0.3;% decay time for E synapses (ms)
%     vars.tau_r_I = 0.5 * mI;% rise time for I synapses (ms)
%     vars.tau_d_I = 2 * mI;% decay time for I synapses (ms)
%
%     vars.mode = 0.5;% the current balance between synaptic weights ( 1:= lateral, 0:= recurrent ), \in [ 0, 1 ].
%     vars.infoStims = [ 4800, 120, 40, 10, 1 ];% recording duration (ms), stimulation duration (ms), stimulation delay (ms), the number of assemblies, the number of stimulations ];
%     vars.inputs = [ 1, 2, 3, 4, 5, 6 ];% in area
%     vars = SNN_stimulation( vars, noiseLevel );
%     [ vars, records ] = SNN_simulation_GPU( vars, 'plasticity_off', 'records_on' );
%
%     s_E = records.s_E;
%     s_E_all{ ct_noiseLevel, 1 } = s_E;
%     s_E = ( transMatU * s_E ) * transMatT;
%     s_E_all_trans{ ct_noiseLevel, 1 } = s_E;
%
%     % figure
%     % imagesc( vars.ext_stim_E ); colorbar
%     % figure
%     % imagesc( records.s_E ); colorbar
%     % figure
%     % imagesc( s_E ); colorbar
%
%
%     Z1 = s_E( 1 : 10, : );
%     Z2 = s_E( 11 : 20, : );
%     Z3 = s_E( 21 : 30, : );
%     Z4 = s_E( 31 : 40, : );
%     Z5 = s_E( 41 : 50, : );
%     Z6 = s_E( 51 : 60, : );
%
%     M1 = mean( sum( Z1, 1 ), 2 );
%     M2 = mean( sum( Z2, 1 ), 2 );
%     M3 = mean( sum( Z3, 1 ), 2 );
%     M4 = mean( sum( Z4, 1 ), 2 );
%     M5 = mean( sum( Z5, 1 ), 2 );
%     M6 = mean( sum( Z6, 1 ), 2 );
%     [ ~, C1, H1 ] = patternMatching( S, Z1 );
%     [ ~, C2, H2 ] = patternMatching( S, Z2 );
%     [ ~, C3, H3 ] = patternMatching( S, Z3 );
%     [ ~, C4, H4 ] = patternMatching( S, Z4 );
%     [ ~, C5, H5 ] = patternMatching( S, Z5 );
%     [ ~, C6, H6 ] = patternMatching( S, Z6 );
%     C1 = C1( H1 > 0 );
%     C2 = C2( H2 > 0 );
%     C3 = C3( H3 > 0 );
%     C4 = C4( H4 > 0 );
%     C5 = C5( H5 > 0 );
%     C6 = C6( H6 > 0 );
%     if length( C1 ) >= 0.5 * length( H1 )
%         t1C = C1( 1 : end - 1 );
%         t2C = C1( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C1 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C1 = false;
%     end
%     if length( C2 ) >= 0.5 * length( H2 )
%         t1C = C2( 1 : end - 1 );
%         t2C = C2( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C2 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C2 = false;
%     end
%     if length( C3 ) >= 0.5 * length( H2 )
%         t1C = C3( 1 : end - 1 );
%         t2C = C3( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C3 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C3 = false;
%     end
%     if length( C4 ) >= 0.5 * length( H2 )
%         t1C = C4( 1 : end - 1 );
%         t2C = C4( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C4 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C4 = false;
%     end
%     if length( C5 ) >= 0.5 * length( H2 )
%         t1C = C5( 1 : end - 1 );
%         t2C = C5( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C5 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C5 = false;
%     end
%     if length( C6 ) >= 0.5 * length( H2 )
%         t1C = C6( 1 : end - 1 );
%         t2C = C6( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C6 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C6 = false;
%     end
%     H1 = mean( H1, 2 );
%     H2 = mean( H2, 2 );
%     H3 = mean( H3, 2 );
%     H4 = mean( H4, 2 );
%     H5 = mean( H5, 2 );
%     H6 = mean( H6, 2 );
%
%     M( ct_noiseLevel, 1 ) = mean( [ M1, M2, M3, M4, M5, M6 ], 2 );
%     C( ct_noiseLevel, 1 ) = all( [ C1, C2, C3, C4, C5, C6 ], 2 );
%     H( ct_noiseLevel, 1 ) = mean( [ H1, H2, H3, H4, H5, H6 ], 2 );
%
%     disp( [ 'NoiseLevel: ', num2str( noiseLevel ) ] )
% end; clear noiseLevel
%
% save( 'SNN_Simul6_2.mat', 'M', 'C', 'H', 's_E_all', 's_E_all_trans', 'pool_E', 'pool_I', '-v7.3' )
%

%% SNN_Simul7_2
%
% transMatU = zeros( 70, 3500 );
% ct_u = 0;
% for u = 1 : 70
%     transMatU( u, ct_u + [ 1 : 50 ] ) = 1;
%     ct_u = ct_u + 50;
% end; clear u ct_u
%
% transMatT = zeros( 9600, 480 );
% ct_t = 0;
% for t = 1 : 480
%     transMatT( ct_t + [ 1 : 20 ], t ) = 1;
%     ct_t = ct_t + 20;
% end; clear t ct_t
%
% n_stim = 10;
% S = eye( n_stim );
%
%
% pool_E = [ 7 : 2 : 25 ];
% pool_I = [ 3 : 3 : 30 ];
%
% M = [];
% C = [];
% H = [];
% s_E_all = {};
% s_E_all_trans = {};
%
% mE = pool_E( 5 );
% mI = pool_I( 7 );
% ct_noiseLevel = 0;
% for noiseLevel = [ 0 : 0.2 : 4 ]
%
%     ct_noiseLevel = ct_noiseLevel + 1;
%
%     load( 'SNN_Simul7_0.mat' )
%
%     vars.tau_r_E = 1 * mE;% rise time for E synapses (ms)
%     vars.tau_d_E = 6 * mE * 0.3;% decay time for E synapses (ms)
%     vars.tau_r_I = 0.5 * mI;% rise time for I synapses (ms)
%     vars.tau_d_I = 2 * mI;% decay time for I synapses (ms)
%
%     vars.mode = 0.5;% the current balance between synaptic weights ( 1:= lateral, 0:= recurrent ), \in [ 0, 1 ].
%     vars.infoStims = [ 4800, 120, 40, 10, 1 ];% recording duration (ms), stimulation duration (ms), stimulation delay (ms), the number of assemblies, the number of stimulations ];
%     vars.inputs = [ 1, 2, 3, 4, 5, 6, 7 ];% in area
%     vars = SNN_stimulation( vars, noiseLevel );
%     [ vars, records ] = SNN_simulation_GPU( vars, 'plasticity_off', 'records_on' );
%
%     s_E = records.s_E;
%     s_E_all{ ct_noiseLevel, 1 } = s_E;
%     s_E = ( transMatU * s_E ) * transMatT;
%     s_E_all_trans{ ct_noiseLevel, 1 } = s_E;
%
%     % figure
%     % imagesc( vars.ext_stim_E ); colorbar
%     % figure
%     % imagesc( records.s_E ); colorbar
%     % figure
%     % imagesc( s_E ); colorbar
%
%
%     Z1 = s_E( 1 : 10, : );
%     Z2 = s_E( 11 : 20, : );
%     Z3 = s_E( 21 : 30, : );
%     Z4 = s_E( 31 : 40, : );
%     Z5 = s_E( 41 : 50, : );
%     Z6 = s_E( 51 : 60, : );
%     Z7 = s_E( 61 : 70, : );
%
%     M1 = mean( sum( Z1, 1 ), 2 );
%     M2 = mean( sum( Z2, 1 ), 2 );
%     M3 = mean( sum( Z3, 1 ), 2 );
%     M4 = mean( sum( Z4, 1 ), 2 );
%     M5 = mean( sum( Z5, 1 ), 2 );
%     M6 = mean( sum( Z6, 1 ), 2 );
%     M7 = mean( sum( Z7, 1 ), 2 );
%     [ ~, C1, H1 ] = patternMatching( S, Z1 );
%     [ ~, C2, H2 ] = patternMatching( S, Z2 );
%     [ ~, C3, H3 ] = patternMatching( S, Z3 );
%     [ ~, C4, H4 ] = patternMatching( S, Z4 );
%     [ ~, C5, H5 ] = patternMatching( S, Z5 );
%     [ ~, C6, H6 ] = patternMatching( S, Z6 );
%     [ ~, C7, H7 ] = patternMatching( S, Z7 );
%     C1 = C1( H1 > 0 );
%     C2 = C2( H2 > 0 );
%     C3 = C3( H3 > 0 );
%     C4 = C4( H4 > 0 );
%     C5 = C5( H5 > 0 );
%     C6 = C6( H6 > 0 );
%     C7 = C7( H7 > 0 );
%     if length( C1 ) >= 0.5 * length( H1 )
%         t1C = C1( 1 : end - 1 );
%         t2C = C1( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C1 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C1 = false;
%     end
%     if length( C2 ) >= 0.5 * length( H2 )
%         t1C = C2( 1 : end - 1 );
%         t2C = C2( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C2 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C2 = false;
%     end
%     if length( C3 ) >= 0.5 * length( H2 )
%         t1C = C3( 1 : end - 1 );
%         t2C = C3( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C3 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C3 = false;
%     end
%     if length( C4 ) >= 0.5 * length( H2 )
%         t1C = C4( 1 : end - 1 );
%         t2C = C4( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C4 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C4 = false;
%     end
%     if length( C5 ) >= 0.5 * length( H2 )
%         t1C = C5( 1 : end - 1 );
%         t2C = C5( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C5 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C5 = false;
%     end
%     if length( C6 ) >= 0.5 * length( H2 )
%         t1C = C6( 1 : end - 1 );
%         t2C = C6( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C6 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C6 = false;
%     end
%     if length( C7 ) >= 0.5 * length( H2 )
%         t1C = C7( 1 : end - 1 );
%         t2C = C7( 2 : end );
%         idxC = [];
%         for n = 1 : n_stim
%             if n < n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n + 1;
%             elseif n == n_stim
%                 idx = t1C == n;
%                 idx = t2C( idx ) == n | t2C( idx ) == n - ( n_stim - 1 );
%             end
%             if ~isempty( idx )
%                 idxC = [ idxC, all( idx, 2 ) ];
%             else
%                 idxC = [ idxC, false ];
%             end
%             clear idx
%         end; clear n
%         C7 = mean( idxC, 2 ) >= 0.5;
%         clear t1C t2C idxC
%     else
%         C7 = false;
%     end
%     H1 = mean( H1, 2 );
%     H2 = mean( H2, 2 );
%     H3 = mean( H3, 2 );
%     H4 = mean( H4, 2 );
%     H5 = mean( H5, 2 );
%     H6 = mean( H6, 2 );
%     H7 = mean( H7, 2 );
%
%     M( ct_noiseLevel, 1 ) = mean( [ M1, M2, M3, M4, M5, M6, M7 ], 2 );
%     C( ct_noiseLevel, 1 ) = all( [ C1, C2, C3, C4, C5, C6, C7 ], 2 );
%     H( ct_noiseLevel, 1 ) = mean( [ H1, H2, H3, H4, H5, H6, H7 ], 2 );
%
%     disp( [ 'NoiseLevel: ', num2str( noiseLevel ) ] )
% end; clear noiseLevel
%
% save( 'SNN_Simul7_2.mat', 'M', 'C', 'H', 's_E_all', 's_E_all_trans', 'pool_E', 'pool_I', '-v7.3' )
%

%% Result, SNN SimulX 2
% 
% noiseLevel = [ 0 : 0.2 : 4 ];
% 
% C_pool = [];
% for k = 2 : 6
%     fname = [ 'SNN_Simul', num2str( k ), '_2.mat' ];
%     load( fname )
%     C_pool = [ C_pool, C ];
% end; clear k
% 
% colors = [ 1, 1, 1; 0.8, 0.8, 0.3 ];
% 
% figure( 'position', [ 100, 100, 300, 150 ] )
% imagesc( noiseLevel, [ 2 : 6 ], transpose( C_pool ) )
% colormap( colors )
% colorbar( 'ticks', [ 0, 1 ], 'ticklabels', { '', '' } )
% 
% % -------------------------------------------------------------------------
% 
% load( 'SNN_Simul3_2.mat' )
% 
% noiseLevel = [ 0 : 0.2 : 4 ];
% 
% % colors = gray;
% % colors = flipud( colors );
% %
% % ct_n = 5;
% % s_E = s_E_all{ ct_n, 1 };
% % figure( 'position', [ 100, 100, 300, 500 ] )
% % imagesc( linspace( 0, 500, 2 * 500 ), [ 1 : 1500 ], s_E( :, 2 * 50 + [ 1 : 2 * 500 ] ) )
% % colormap( colors )
% % disp( [ 'noiseLevel = ', num2str( noiseLevel( ct_n ) ) ] )
% %
% % ct_n = 7;
% % s_E = s_E_all{ ct_n, 1 };
% % figure( 'position', [ 100, 100, 300, 500 ] )
% % imagesc( linspace( 0, 500, 2 * 500 ), [ 1 : 1500 ], s_E( :, 2 * 50 + [ 1 : 2 * 500 ] ) )
% % colormap( colors )
% % disp( [ 'noiseLevel = ', num2str( noiseLevel( ct_n ) ) ] )
% 
% ct_n = 1;
% s_E = s_E_all{ ct_n, 1 };
% figure( 'position', [ 100, 100, 300, 500 ] )
% hold on
% for n = 1 : 500
%     pts = find( s_E( n, 2 * 50 + [ 1 : 2 * 500 ] ) == 1 );
%     if ~isempty( pts )
%         plot( pts, n, '.k' )
%     end
% end; clear n pts
% set( gca, 'xlim', [ 0, 1001 ], 'ylim', [ 0, 501 ] )
% set( gca, 'xtick', linspace( 1, 1000, 6 ), 'xticklabel', { '0', '100', '200', '300', '400', '500' } )
% disp( [ 'noiseLevel = ', num2str( noiseLevel( ct_n ) ) ] )
% 
% ct_n = 5;
% s_E = s_E_all{ ct_n, 1 };
% figure( 'position', [ 100, 100, 300, 500 ] )
% hold on
% for n = 1 : 500
%     pts = find( s_E( n, 2 * 50 + [ 1 : 2 * 500 ] ) == 1 );
%     if ~isempty( pts )
%         plot( pts, n, '.k' )
%     end
% end; clear n pts
% set( gca, 'xlim', [ 0, 1001 ], 'ylim', [ 0, 501 ] )
% set( gca, 'xtick', linspace( 1, 1000, 6 ), 'xticklabel', { '0', '100', '200', '300', '400', '500' } )
% disp( [ 'noiseLevel = ', num2str( noiseLevel( ct_n ) ) ] )
% 
% ct_n = 6;
% s_E = s_E_all{ ct_n, 1 };
% figure( 'position', [ 100, 100, 300, 500 ] )
% hold on
% for n = 1 : 500
%     pts = find( s_E( n, 2 * 50 + [ 1 : 2 * 500 ] ) == 1 );
%     if ~isempty( pts )
%         plot( pts, n, '.k' )
%     end
% end; clear n pts
% set( gca, 'xlim', [ 0, 1001 ], 'ylim', [ 0, 501 ] )
% set( gca, 'xtick', linspace( 1, 1000, 6 ), 'xticklabel', { '0', '100', '200', '300', '400', '500' } )
% disp( [ 'noiseLevel = ', num2str( noiseLevel( ct_n ) ) ] )
% 
