function [ W1, W2, W3, W4, W5, W6 ] = learningHebbian6( S1, S2, S3, S4, S5, S6 )

W1 = zeros( size( S1, 1 ), size( S1, 1 ) );
W2 = zeros( size( S2, 1 ), size( S2, 1 ) );
W3 = zeros( size( S3, 1 ), size( S3, 1 ) );
W4 = zeros( size( S4, 1 ), size( S4, 1 ) );
W5 = zeros( size( S5, 1 ), size( S5, 1 ) );
W6 = zeros( size( S6, 1 ), size( S6, 1 ) );
for t = 1 : size( S1, 2 )
    W1( S1( :, t ) == 1, S2( :, t ) == 1) = W1( S1( :, t ) == 1, S2( :, t ) == 1) + 1;
    W2( S2( :, t ) == 1, S3( :, t ) == 1) = W2( S2( :, t ) == 1, S3( :, t ) == 1) + 1;
    W3( S3( :, t ) == 1, S4( :, t ) == 1) = W3( S3( :, t ) == 1, S4( :, t ) == 1) + 1;
    W4( S4( :, t ) == 1, S5( :, t ) == 1) = W4( S4( :, t ) == 1, S5( :, t ) == 1) + 1;
    W5( S5( :, t ) == 1, S6( :, t ) == 1) = W5( S5( :, t ) == 1, S6( :, t ) == 1) + 1;
    W6( S6( :, t ) == 1, S1( :, t ) == 1) = W6( S6( :, t ) == 1, S1( :, t ) == 1) + 1;
end
W1 = bsxfun( @rdivide, W1, sum( W1, 2 ) );
W2 = bsxfun( @rdivide, W2, sum( W2, 2 ) );
W3 = bsxfun( @rdivide, W3, sum( W3, 2 ) );
W4 = bsxfun( @rdivide, W4, sum( W4, 2 ) );
W5 = bsxfun( @rdivide, W5, sum( W5, 2 ) );
W6 = bsxfun( @rdivide, W6, sum( W6, 2 ) );
