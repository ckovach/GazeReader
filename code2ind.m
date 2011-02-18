function ind = code2ind(h,code)


regData = getappdata(h,'regData');

codes = [regData.regressors.code];

[q,ind] = ismember(code,codes);
