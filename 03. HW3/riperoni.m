%non-unique reduction, don't try to make sense out of this
clc, clear, close all

% A = [ 13 13 24 25 30 ; 2 3 2 3 5 ];
% B = [ 2 7 1 3 10 ];

matches = [1 6 7 8 9 11 13 15 16 18 19 24; 14 32 24 25 6 8 32 14 32 10 8 32];
scores = [ 64356 54119 66032 76157 99879 150982 64073 62596 75491 116415 207530 131190];

function nonUniques(matches,scores)
Aa = matches(1,:);
Ab = matches(2,:);
Bc = scores;

n = length(Bc);3
while n ~= 0
    nn = 0;
    for ii = 1:n-1
        if Bc(ii) > Bc(ii+1)
            [Bc(ii+1),Bc(ii)] = deal(Bc(ii), Bc(ii+1));
            [Aa(ii+1),Aa(ii)] = deal(Aa(ii), Aa(ii+1));
            [Ab(ii+1),Ab(ii)] = deal(Ab(ii), Ab(ii+1));
            nn = ii;
        end
    end
    n = nn;
end

A2 = [Aa;Ab];
B2 = Bc;

A3 = [];
B3 = [];
k = 1;
for j = 1:length(B2)           %Reduce to best unique matches
    P = A2(1,j);
    if(isempty(find(A3 == P)))
        Q = A2(2,j);
        idxOfUniqueMinScore = find(B2 == min(B2(round((find(A2 == Q).')/2))));
        idxOfUniqueMinMatch = A2(2,idxOfUniqueMinScore(1));
        if(isempty(find(A3 == Q)))
            A3(1,k) = A2(1,j);
            A3(2,k) = A2(2,idxOfUniqueMinScore);
            B3(k) = B2(idxOfUniqueMinScore(1));
            k = k + 1;
        end
        
    end
    
end
end