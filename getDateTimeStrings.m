function [dateStr, timeStr] = getDateTimeStrings(my_date,my_clock)

dateStr = regexprep(my_date,'-','');
timeStr1_entries = cellfun(@(x) num2str(x),num2cell(fix(my_clock(4:6))),'UniformOutput',false);
for i=1:size(timeStr1_entries,2)
    if numel(timeStr1_entries{i}) == 1
        timeStr1_entries{i} = cell2mat(['0',timeStr1_entries(i)]);
    end
end
timeStr = cell2mat(timeStr1_entries);

end