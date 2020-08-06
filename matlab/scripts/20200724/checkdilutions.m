% Check dilution series
dildata=dir([msdir,'20200724/CDIV-*.mzXML']);
initmap=struct('time',[0,197;3300,2078],'mz',[ 133,133-17e-4;500,500-79e-4]);
ref=31;
if ~exist('dilms','var')
  dilms={};
end
if ~exist('dilmaps','var')
  dilmaps=[];
end
for i=1:length(dildata)
  if i> length(dilms)
    path=[dildata(i).folder,'/',dildata(i).name];
    dilms{i}=MassSpec(path);
    dilms{i}.setLoad(2000*1e-15);
    % Prune out some data
    dilms{i}.filter([200,2900],[130,530]);  % NOTE: this is using the localtimes to filter
  end
  if i>length(dilmaps)
    dilmaps=[dilmaps,compounds.computeMap(dilms{i},'initmap',initmap)];
  end
  compounds.plotComposition(dilms{i},'mztol',.01,'map',dilmaps(i),'ref',mzdata{ref},'refmap',maps(ref));
end

mz=compounds.mass+compounds.ADDUCTS(1).mass;
mz=mz(:);
for i=[1,3,4]
  dilms{2}.comparepeaks(dilms{i},'mztol',.01,'minic',5000,'mz',mz);
end
full=mzdata{31};
full.comparepeaks(dilms{2},'mztol',.01,'minic',10000,'mz',mz);
oldcdiv=mzdata{60};
full.comparepeaks(oldcdiv,'mztol',.01,'minic',10000,'mz',mz);
f8630=mzdata{44};  
full.comparepeaks(f8630,'mztol',.01,'minic',10000,'mz',mz);
