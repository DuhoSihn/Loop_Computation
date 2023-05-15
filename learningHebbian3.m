function [ W1, W2, W3 ] = learningHebbian3( S1, S2, S3 )

W1 = zeros( size( S1, 1 ), size( S1, 1 ) );
W2 = zeros( size( S2, 1 ), size( S2, 1 ) );
W3 = zeros( size( S3, 1 ), size( S3, 1 ) );
for t = 1 : size( S1, 2 )
    W1( S1( :, t ) == 1, S2( :, t ) == 1) = W1( S1( :, t ) == 1, S2( :, t ) == 1) + 1;
    W2( S2( :, t ) == 1, S3( :, t ) == 1) = W2( S2( :, t ) == 1, S3( :, t ) == 1) + 1;
    W3( S3( :, t ) == 1, S1( :, t ) == 1) = W3( S3( :, t ) == 1, S1( :, t ) == 1) + 1;
end
W1 = bsxfun( @rdivide, W1, sum( W1, 2 ) );
W2 = bsxfun( @rdivide, W2, sum( W2, 2 ) );
W3 = bsxfun( @rdivide, W3, sum( W3, 2 ) );
