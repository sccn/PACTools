function studyflag = isstudy(struct)
studyflag = 0;
if isstruct(struct)
    Allfields = fieldnames(struct);
    studyflag = any(find(~cellfun(@isempty,strfind(Allfields,'cluster' )))) && any(find(~cellfun(@isempty,strfind(Allfields,'datasetinfo' ))));
end