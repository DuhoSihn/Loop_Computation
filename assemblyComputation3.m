function [ Z1, Z2, Z3, W1, W2, W3 ] = assemblyComputation3( S1, S2, S3, W1, W2, W3, threshold, period_initiation, period_active, period_refractory, learning_rate )
% [ Z1, Z2, Z3 ] = assemblyComputation3( S1, S2, S3, W1, W2, W3, threshold, period_initiation, period_active, period_refractory )

if size( S1, 2 ) ~= size( S2, 2 ), error( 'Duration of S1 and S2 must be same!' ), end
if size( S2, 2 ) ~= size( S3, 2 ), error( 'Duration of S2 and S3 must be same!' ), end

Z1 = NaN( size( S1 ) );
Z2 = NaN( size( S2 ) );
Z3 = NaN( size( S3 ) );

count1 = NaN( size( S1 ) );
count2 = NaN( size( S2 ) );
count3 = NaN( size( S3 ) );

for t = 1 : size( S1, 2 )
    
    if any( isnan( S1( :, t ) ), 1 )
        X1 = Z1( :, t - 1 );
    else
        X1 = S1( :, t );
    end
    if any( isnan( S2( :, t ) ), 1 )
        X2 = Z2( :, t - 1 );
    else
        X2 = S2( :, t );
    end
    if any( isnan( S3( :, t ) ), 1 )
        X3 = Z3( :, t - 1 );
    else
        X3 = S3( :, t );
    end
    
    Y1 = W1 * X2;
    Y2 = W2 * X3;
    Y3 = W3 * X1;
    Y1 = bsxfun( @ge, Y1, threshold );
    Y2 = bsxfun( @ge, Y2, threshold );
    Y3 = bsxfun( @ge, Y3, threshold );
    
    if t == 1
        
        count1( :, t ) = 0;
        count2( :, t ) = 0;
        count3( :, t ) = 0;
        
        idx1 = count1( :, t ) == 0 & Y1 == 1;
        idx2 = count2( :, t ) == 0 & Y2 == 1;
        idx3 = count3( :, t ) == 0 & Y3 == 1;
        count1( idx1, t ) = 1;
        count2( idx2, t ) = 1;
        count3( idx3, t ) = 1;
        
    else
        
        count1( :, t ) = 0;
        count2( :, t ) = 0;
        count3( :, t ) = 0;
        
        idx1 = count1( :, t - 1 ) >= 1 & count1( :, t - 1 ) < period_initiation + period_active + period_refractory;
        idx2 = count2( :, t - 1 ) >= 1 & count2( :, t - 1 ) < period_initiation + period_active + period_refractory;
        idx3 = count3( :, t - 1 ) >= 1 & count3( :, t - 1 ) < period_initiation + period_active + period_refractory;
        count1( idx1, t ) = count1( idx1, t - 1 ) + 1;
        count2( idx2, t ) = count2( idx2, t - 1 ) + 1;
        count3( idx3, t ) = count3( idx3, t - 1 ) + 1;
        
        idx1 = count1( :, t - 1 ) == period_initiation + period_active + period_refractory;
        idx2 = count2( :, t - 1 ) == period_initiation + period_active + period_refractory;
        idx3 = count3( :, t - 1 ) == period_initiation + period_active + period_refractory;
        count1( idx1, t ) = 0;
        count2( idx2, t ) = 0;
        count3( idx3, t ) = 0;
        
        idx1 = count1( :, t ) == 0 & Y1 == 1;
        idx2 = count2( :, t ) == 0 & Y2 == 1;
        idx3 = count3( :, t ) == 0 & Y3 == 1;
        count1( idx1, t ) = 1;
        count2( idx2, t ) = 1;
        count3( idx3, t ) = 1;
        
    end
    
    Z1( :, t ) = 0;
    Z2( :, t ) = 0;
    Z3( :, t ) = 0;
    idx1 = count1( :, t ) > period_initiation & count1( :, t ) <= period_initiation + period_active;
    idx2 = count2( :, t ) > period_initiation & count2( :, t ) <= period_initiation + period_active;
    idx3 = count3( :, t ) > period_initiation & count3( :, t ) <= period_initiation + period_active;
    Z1( idx1, t ) = 1;
    Z2( idx2, t ) = 1;
    Z3( idx3, t ) = 1;
    
    if nargin == 9
        W1( X1 == 1, X2 == 1) = ( 1 + learning_rate ) * W1( X1 == 1, X2 == 1);
        W2( X2 == 1, X3 == 1) = ( 1 + learning_rate ) * W2( X2 == 1, X3 == 1);
        W3( X3 == 1, X1 == 1) = ( 1 + learning_rate ) * W3( X3 == 1, X1 == 1);
        W1 = bsxfun( @rdivide, W1, sum( W1, 2 ) );
        W2 = bsxfun( @rdivide, W2, sum( W2, 2 ) );
        W3 = bsxfun( @rdivide, W3, sum( W3, 2 ) );
    end
    
end
