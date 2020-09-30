function newStruct = sameFields(oldStruct,thisFields)
%Make sure the struct has the same fields
newStruct = oldStruct;
oldFields = fields(oldStruct);
for ii=1:length(oldFields)
    if ~ismember(oldFields{ii}, thisFields)
        newStruct = rmfield(newStruct, oldFields{ii});
    end
end

end

