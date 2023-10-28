function [ W1, W2 ] = learningHebbian2( S1, S2 )

W1 = zeros( size( S1, 1 ), size( S1, 1 ) );
W2 = zeros( size( S2, 1 ), size( S2, 1 ) );
for t = 1 : size( S1, 2 )
    W1( S1( :, t ) == 1, S2( :, t ) == 1) = W1( S1( :, t ) == 1, S2( :, t ) == 1) + 1;
    W2( S2( :, t ) == 1, S1( :, t ) == 1) = W2( S2( :, t ) == 1, S1( :, t ) == 1) + 1;
end
W1 = bsxfun( @rdivide, W1, sum( W1, 2 ) );
W2 = bsxfun( @rdivide, W2, sum( W2, 2 ) );
