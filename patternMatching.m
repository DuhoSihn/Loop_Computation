function [ D, C, H ] = patternMatching( S, Z )

% distance from S to Z
D = permute( sqrt( sum( bsxfun( @minus, Z, permute( S, [ 1, 3, 2 ] ) ) .^ 2, 1 ) ), [ 3, 2, 1 ] );

% class
[ minD, C ] = min( D, [], 1 );

% certainty
meanD = ( sum( D, 1 ) - minD ) / ( size( S, 2 ) - 1 );
H = meanD - minD;
