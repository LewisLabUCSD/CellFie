function longm = wide2long(ds,rownames,colnames)
	% dstmp=mat2dataset(ds)
	% dstmp.row = rownames(:,2)
	% dsNew = stack(dstmp,dstmp.Properties.VarNames(1:end-1),...
 %            'newDataVarName','Scores')
	%ds.Properties.VarNames = colnames
	longm = {};
	longm.sample = [];
	longm.task = {};
	longm.value = [];
	count=1;
	for i=1:size(ds,2)
		for j=1:size(ds,1)
			longm.sample(count) = colnames(i);
			longm.task{count} = rownames{j};
			longm.value(count)= ds(j,i);
			count=count+1;
        end
    end
end
    