%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = blocking(A,B,delta)
%BLOCKING  Produce blocking pattern for block Parlett recurrence in FUNM.
%   M = BLOCKING(A, DELTA, SHOWPLOT) accepts an upper triangular matrix
%   A and produces a blocking pattern, specified by the vector M,
%   for the block Parlett recurrence.
%   M(i) is the index of the block into which A(i,i) should be placed,
%   for i=1:LENGTH(A).
%   DELTA is a gap parameter (default 0.1) used to determine the blocking.

%   For A coming from a real matrix it should be posible to take
%   advantage of the symmetry about the real axis.  This code does not.

%   References:
%   P. I. Davies and N. J. Higham, A Schur-Parlett algorithm for computing
%      matrix functions. SIAM J. Matrix Anal. Appl., 25(2):464-485, 2003.
%   N. J. Higham, Functions of Matrices: Theory and Computation,
%      Society for Industrial and Applied Mathematics, Philadelphia, PA,
%      USA, 2008.
%
%   Nicholas J. Higham
%   Copyright 1984-2017 The MathWorks, Inc. 


a = diag(A); n = length(a);
b = diag(B);
m = zeros(1,n); maxM = 0;

if nargin < 3 || isempty(delta), delta = 1e-10; end

for i = 1:n

    if m(i) == 0
        m(i) = maxM + 1; % If a(i) hasn`t been assigned to a set
        maxM = maxM + 1; % then make a new set and assign a(i) to it.
    end

    for j = i+1:n
        if m(i) ~= m(j)    % If a(i) and a(j) are not in same set.
            if abs(a(i)*b(j)-b(i)*a(j)) <= delta*((abs(a(i))^2+abs(b(i))^2)^(1/2)*(abs(a(j))^2+abs(b(j))^2)^(1/2))          %%%AMG 9/11/19

                if m(j) == 0
                    m(j) = m(i); % If a(j)/b(j) hasn`t been assigned to a %%%AMG 9/11/19
                                 % set, assign it to the same set as a(i)/b(i).%%%AMG 9/11/19
                else
                    p = max(m(i),m(j)); q = min(m(i),m(j));
                    m(m==p) = q; % If a(j)/b(j) has been assigned to a set %%%AMG 9/11/19
                                 % place all the elements in the set
                                 % containing a(j)/b(j) into the set %%%AMG 9/11/19
                                 % containing a(i)/b(j) (or vice versa).  %%%AMG 9/11/19
                    m(m>p) = m(m>p) -1;
                    maxM = maxM - 1;
                                 % Tidying up. As we have deleted set
                                 % p we reduce the index of the sets
                                 % > p by 1.
                end
            end
        end
    end
end
end