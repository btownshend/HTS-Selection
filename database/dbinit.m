% Populate compounds database from HTBC 5120 set

datadir='../../data/';
matdir=[datadir,'matfiles/'];
resultsdir='../../results/';

if ~exist('s7sdf','var')
  load([matdir,'s7sdf.mat']);
end

if ~exist('s7vecs','var')
  s7vecs=load([matdir,'s7vecs.mat']);
  s7vecs=s7vecs.s7vecs;
end

if ~exist('s7inchi','var')
  % Mapping from STF numbers to InChIKeys from David 2/2021
  s7inchi=readtable([datadir,'HTBCFiles/SmolkeCDIV.csv']);
end

% Create CDIV and CDIQ naming domains
if mysql('status')~=0
  dbhost='35.203.151.202';
  dbuser='hts';
  dbpassword='driver';
  mysql('open',dbhost,dbuser,dbpassword);
end
mysql('use compounds');
cdivdomain=mysql('select domain from domains where name=''CDIV''');
if isempty(cdivdomain)
  mysql('insert into domains(name) values(''CDIV'')');
  cdivdomain=mysql('select domain from domains where name=''CDIV''');
end
cdiqdomain=mysql('select domain from domains where name=''CDIQ''');
if isempty(cdiqdomain)
  mysql('insert into domains(name) values(''CDIQ'')');
  cdiqdomain=mysql('select domain from domains where name=''CDIQ''');
end
stfdomain=mysql('select domain from domains where name=''STF''');
if isempty(stfdomain)
  mysql('insert into domains(name) values(''STF'')');
  stfdomain=mysql('select domain from domains where name=''STF''');
end
vendordomain=mysql('select domain from domains where name=''ChemDiv''');
if isempty(vendordomain)
  mysql('insert into domains(name) values(''ChemDiv'')');
  vendordomain=mysql('select domain from domains where name=''ChemDiv''');
end

mysql('delete from names');
mysql('delete from contents');
mysql('delete from mixtures');
mysql('delete from compounds');
mysql('ALTER TABLE compounds AUTO_INCREMENT=1');
mysql('ALTER TABLE mixtures AUTO_INCREMENT=1');
mysql('insert into compounds(name,inChIKey,formula,monoisotopicMass) values(''Water'',''XLYOFNOQVPJJNP-UHFFFAOYSA-N'',''H2O'',18.010564684)');
water=mysql('select LAST_INSERT_ID()');
mysql('insert into compounds(name,inChIKey,formula,monoisotopicMass) values(''DMSO'',''IAZDPXIOMUYVGZ-WFGJKAKNSA-N'',''C2H6OS'',84.05159646)');
dmso=mysql('select LAST_INSERT_ID()');

aliases={};
for i=1:length(s7sdf.sdf)
  s=s7sdf.sdf(i);
  formula=s7sdf.getformula(i);
  stf=str2num(s.compound_Corp_Reg_Number(5:end));
  cdiv=sprintf('%d%s',str2num(s.BATCH_PLATE(5:end)),s.BATCH_WELL);
  cdiq=sprintf('%d%s',str2num(s.Plate384(5:end)),s.Well384);
  ind1=strcmp([s.BATCH_PLATE,'-',s.BATCH_WELL],s7vecs.cdivnames);
  assert(sum(ind1)==1);
  assert(strcmp([s.Plate384,'-',s.Well384],s7vecs.cdiqnames{ind1}));
  ind2=strcmp(s7inchi.ID,sprintf('STF-%08d',stf));
  assert(sum(ind2)==1);
  inchikey=s7inchi.InChIKey{ind2};
  vendor=s7inchi.SUPPLIER_REF1{ind2};
  mmass=s.MonoisotopicMass;
  mysql(sprintf('INSERT INTO compounds(name,inChIKey,monoisotopicMass,formula) values(''%s'',''%s'',%f,''%s'')',cdiv,inchikey,mmass,formula));
  cid=mysql('select LAST_INSERT_ID()');
  fprintf('CID=%d: %d %d STF-%d CDIV-%s CDIQ-%s %s %s %f\n', cid, i, find(ind1), stf, cdiv, cdiq, vendor, inchikey, mmass);
  aliases{end+1}=sprintf('(%d,%d,''%s'')',cid,cdivdomain,cdiv);
  aliases{end+1}=sprintf('(%d,%d,''%s'')',cid,cdiqdomain,cdiq);
  aliases{end+1}=sprintf('(%d,%d,''%d'')',cid,stfdomain,stf);
  aliases{end+1}=sprintf('(%d,%d,''%s'')',cid,vendordomain,vendor);
end

% Add all the aliases
mysql(sprintf('INSERT INTO names(compound,domain,name) values %s',strjoin(aliases,',')));

% Add all the CDIQ plates
cmd=['INSERT INTO mixtures(name,solvent) ',...
'SELECT CONCAT(''CDIQ-'',n.name),sv.compound ',...
'FROM names n, compounds sv ',...
'WHERE n.domain=(SELECT domain FROM domains d WHERE d.name=''CDIQ'') ',...
'AND sv.name=''DMSO'' ']
mysql(cmd);

cmd=['INSERT INTO contents(mixture,compound,concentration) ',...
'SELECT m.mixture, n.compound, 100e-6 ',...
'FROM mixtures m, names n ',...
'WHERE m.name=CONCAT(''CDIQ-'',n.name) ',...
'AND n.domain=(SELECT domain FROM domains d WHERE d.name=''CDIQ'')']
mysql(cmd);

% Add the R320 vectors
cmd=['INSERT INTO mixtures(name,solvent) ',...
'SELECT DISTINCT CONCAT(''R320-'',LEFT(RIGHT(n.name,3),1)),sv.compound ',...
'FROM names n, compounds sv ',...
'WHERE n.domain=(SELECT domain FROM domains d WHERE d.name=''CDIQ'') ',...
'AND sv.name=''DMSO'' ']
mysql(cmd);

cmd=['INSERT INTO contents(mixture,compound,concentration) ',...
'SELECT m.mixture,n.compound, 5e-3/320 ',...
'FROM mixtures m, names n ',...
'WHERE n.domain=(SELECT domain FROM domains d WHERE d.name=''CDIQ'') ',...
'AND CONCAT(''R320-'',LEFT(RIGHT(n.name,3),1))=m.name ']
mysql(cmd);

% Add V2560A/B
cmd=['INSERT INTO mixtures(name,solvent) ',...
'SELECT ''V2560A'',sv.compound ',...
'FROM compounds sv ',...
'WHERE sv.name=''DMSO'' ']
mysql(cmd);

cmd=['INSERT INTO contents(mixture,compound,concentration) ',...
'SELECT m.mixture,n.compound, 5e-3/2560 ',...
'FROM mixtures m, names n ',...
'WHERE n.domain=(SELECT domain FROM domains d WHERE d.name=''CDIQ'') ',...
'AND LEFT(RIGHT(n.name,3),1)<=''H'' ',...
'AND m.name=''V2560A'' ']
mysql(cmd);

cmd=['INSERT INTO mixtures(name,solvent) ',...
'SELECT ''V2560B'',sv.compound ',...
'FROM compounds sv ',...
'WHERE sv.name=''DMSO'' ']
mysql(cmd);

cmd=['INSERT INTO contents(mixture,compound,concentration) ',...
'SELECT m.mixture,n.compound, 5e-3/2560 ',...
'FROM mixtures m, names n ',...
'WHERE n.domain=(SELECT domain FROM domains d WHERE d.name=''CDIQ'') ',...
'AND LEFT(RIGHT(n.name,3),1)>''H'' ',...
'AND m.name=''V2560B'' ']
mysql(cmd);

% Add V5120
cmd=['INSERT INTO mixtures(name,solvent) ',...
'SELECT ''V5120'',sv.compound ',...
'FROM compounds sv ',...
'WHERE sv.name=''DMSO'' ']
mysql(cmd);

cmd=['INSERT INTO contents(mixture,compound,concentration) ',...
'SELECT m.mixture,n.compound, 5e-3/5120 ',...
'FROM mixtures m, names n ',...
'WHERE n.domain=(SELECT domain FROM domains d WHERE d.name=''CDIQ'') ',...
'AND m.name=''V5120'' ']
mysql(cmd);

% Add the V64 vectors
for i=1:length(s7vecs.v64names)
  name=s7vecs.v64names{i};
  contains=s7vecs.cdiqnames(s7vecs.v64(i,:));
  mysql(sprintf('INSERT INTO mixtures(name,solvent) SELECT ''%s'',compound FROM compounds c WHERE c.name=''DMSO''',name));
  %mix=mysql('select LAST_INSERT_ID()');
  contents=sprintf('INSERT INTO contents(mixture,compound,concentration)\nSELECT LAST_INSERT_ID(),compound,%g FROM names n\n',5e-3/64);
  contents=[contents,sprintf('WHERE domain=%d\nAND n.name IN ',cdiqdomain)];
  nlist=cellfun(@(z) sprintf('''%d%s''',str2num(z(5:8)),z(10:end)),contains,'Unif',false);
  contents=[contents,'(',strjoin(nlist,','),')'];
  nins=mysql(contents)
end

% Add the V256 vectors
for i=1:length(s7vecs.v256names)
  name=s7vecs.v256names{i};
  contains=s7vecs.cdiqnames(s7vecs.v256(i,:));
  mysql(sprintf('INSERT INTO mixtures(name,solvent) SELECT ''%s'',compound FROM compounds c WHERE c.name=''DMSO''',name));
  %mix=mysql('select LAST_INSERT_ID()');
  contents=sprintf('INSERT INTO contents(mixture,compound,concentration)\nSELECT LAST_INSERT_ID(),compound,%g FROM names n\n',5e-3/256);
  contents=[contents,sprintf('WHERE domain=%d\nAND n.name IN ',cdiqdomain)];
  nlist=cellfun(@(z) sprintf('''%d%s''',str2num(z(5:8)),z(10:end)),contains,'Unif',false);
  contents=[contents,'(',strjoin(nlist,','),')'];
  nins=mysql(contents)
end

return;

% Subsequently repaired some formula/masses based on PubChem records:
update compounds set formula='C18H21N2O3S2+',monoisotopicMass=377.09936 where compound=465;  46G04
update compounds set formula='C11H9NO2S',monoisotopicMass=219.0354  where compound=862; 87F11
update compounds set formula='C15H16NO3S2+',monoisotopicMass=322.057161  where compound=1448; 167A07
update compounds set formula='C11H16NO3S3+',monoisotopicMass=306.029232   where compound=1497; 167F06
update compounds set formula='C11H16NO3S3+',monoisotopicMass=306.029232   where compound=1507; 167G06
update compounds set formula='C11H15ClNO3S2+',monoisotopicMass=308.018188   where compound=1517;  167H06


update mixtures
set name=concat(left(name,4),right(left(name,8),3),'-',right(name,3))  
where name like'CDIQ-%'
and length(name)=11;

update mixtures
set name=concat(left(name,4),'0',right(left(name,7),2),'-',right(name,3))  
where name like'CDIQ-%'
and length(name)=10;

update mixtures
set name=concat(left(name,4),'00',right(left(name,6),1),'-',right(name,3))  
where name like'CDIQ-%'
and length(name)=9;
