function ind = code2ind(h,code)


regData = getappdata(h,'regData');

currdat = getappdata(h,'CurrentDataSet');
codes = [regData(currdat).regressors.code];

[q,ind] = ismember(code,codes);
