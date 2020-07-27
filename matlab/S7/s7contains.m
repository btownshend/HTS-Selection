% Figure out what a sample contains based on its name (for S7 5120 set)
function c=s7contains(compounds, vecs, name)
  c=false(1,5120);
  if strncmp(name,'DMSO',4)
    ;   % Nothing
  elseif strcmp(name,'V5120')
    c(:)=true;
  else
    % Line up vectors with compounds 
    vplate96=vecs.csvdata.srcplate96;
    vwell96=vecs.csvdata.srcwell96;
    ord=cellfun(@(z) find(strcmp(z.BATCH_PLATE,vplate96) & strcmp(z.BATCH_WELL,vwell96)), compounds.sdf);
    assert(length(vplate96)==length(ord));
    
    n=strsplit(name,'-');
    assert(length(n)==2);
    row=n{2}(1)-'A'+1;
    if length(n{2})>1
      col=str2num(n{2}(2:end));
    else
      col=nan;
    end
    if strcmp(n{1},'V256A')
      index=(col-1)*8+row;
      c=vecs.v256(index,ord);
      assert(sum(c)==256);
    elseif strcmp(n{1},'V256B')
      index=(col-1)*8+row+96;
      c=vecs.v256(index,ord);
      assert(sum(c)==256);
    elseif strcmp(n{1},'V64A')
      index=(col-1)*16+row;
      c=vecs.v64(index,ord);
      assert(sum(c)==64);
    elseif strcmp(n{1},'V64B')
      index=(col-1)*16+row+384;
      c=vecs.v64(index,ord);
      assert(sum(c)==64);
    elseif strcmp(n{1},'PR80')
      assert(mod(col,3)~=0);
      group=floor((col-1)/3)+1;
      plates=arrayfun(@(z) sprintf('CDIQ%04d',z),(group-1)*160+5:40:group*160+5-1,'Unif',false);
      index=(col-1-(group-1)*3)*8+row;
      c=cellfun(@(z) z.Well384(1)-'A'+1==index & ismember(z.Plate384,plates),compounds.sdf);
      assert(sum(c)==80);
    elseif strcmp(n{1},'R320')
      index=row;
      c=cellfun(@(z) z.Well384(1)-'A'+1==row,compounds.sdf);
      assert(sum(c)==320);
    elseif strncmp(n{1},'CDIQ',4)
      plate=str2num(n{1}(5:end));
      plate384=sprintf('CDIQ%04d',plate);
      well384=sprintf('%c%02d',row+'A'-1,col);
      c=cellfun(@(z) strcmp(z.Plate384,plate384) & strcmp(z.Well384,well384),compounds.sdf);
      assert(sum(c)==1);
    elseif strncmp(n{1},'CDIV',4)
      plate=str2num(n{1}(5:end));
      name=sprintf('%d%c%02d',plate,row+'A'-1,col);
      c=strcmp(compounds.names,name);
      assert(sum(c)==1);
    else
      error('Unable to parse %s',name);
    end
  end
end

