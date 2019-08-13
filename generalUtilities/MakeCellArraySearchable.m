function fixed_cellarray = MakeCellArraySearchable(original_cellarray,replace_str)
% Simple Script to search Cell Array of Strings for numerical entries which
% are converted to variable names (with genvarname) and empty entries which
% are replaced with the string provided in the input "replace_str".
% Default value is "Empty_Cell"

if nargin < 2
    replace_str = 'Empty_Cell';
end

for i=1:size(original_cellarray,1)
    if ~isstr(original_cellarray{i})
        try
            original_cellarray{i} = num2str(original_cellarray{i});
        catch UnknownVarType
            if isempty(original_cellarray{i})
                original_cellarray{i} = replace_str;
            else
                original_cellarray{i} = genvarname(original_cellarray{i});
            end
        end
    end
end

fixed_cellarray = original_cellarray;