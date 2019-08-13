function c = combntns(choicevec,choose)
%COMBNTNS  All possible combinations of a set of values
%
%  COMBNTNS will be removed in a future release. Use NCHOOSEK instead. 
%
%  c = COMBNTNS(choicevec,choose) returns all combinations of the
%  values of the input choice vector.  The size of the combinations
%  are given by the second input.  For example, if choicevec
%  is [1 2 3 4 5], and choose is 2, the output is a matrix
%  containing all distinct pairs of the choicevec set.
%  The output matrix has "choose" columns and the combinatorial
%  "length(choicevec)-choose-'choose'" rows.  The function does not
%  account for repeated values, treating each entry as distinct.
%  As in all combinatorial counting, an entry is not paired with
%  itself, and changed order does not constitute a new pairing.
%  This function is recursive.
%
%  See also NCHOOSEK.

% Copyright 1996-2013 The MathWorks, Inc.

warning(message('map:removing:combntns', 'COMBNTNS'))
c = nchoosek(choicevec, choose);
