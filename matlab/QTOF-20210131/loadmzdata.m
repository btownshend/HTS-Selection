% Setup compounds with elution times for each compound
% Do this by distinguishing isomers by checking row, col, and plate totals separately -- as long as a compound is a unique isomer in at least one of these, then it should be distinguishable

datadir='../../data/';
matdir=[datadir,'matfiles/'];

if ~exist('s7vecs','var')
  load([matdir,'s7vecs.mat']);
end

if ~exist('qsetup','var')
  qsetup=Compounds(.0005,40/60);
  % qsetup.ADDUCTS=qsetup.ADDUCTS(1);   % Only M+H for now
  %qsetup.addCompoundsFromSDF(s7sdf);
  qsetup.addCompoundsFromDB();
end


msdir2=[datadir,'MassSpec/20201002'];
msdir3=[datadir,'MassSpec/2021.01.31'];
data={ % Well, mixture, name, filename, conc, inject
       'A1','DMSO','',[msdir2,'/A1.mzXML'],0,10e-6,
       'B1','V5120','',[msdir2,'/B1.mzXML'],25e-9,10e-6,
       'C1','V5120','',[msdir2,'/C1.mzXML'],100e-9,10e-6,
       'D1','V2560A','',[msdir2,'/D1.mzXML'],100e-9,10e-6,
       'E1','V2560B','',[msdir2,'/E1.mzXML'],100e-9,10e-6,
       'F1','V2560A','V2560A-sn',[msdir2,'/F1.mzXML'],100e-9,10e-6,
       'G1','V2560B','V2560B-sn',[msdir2,'/G1.mzXML'],100e-9,10e-6,
       'H1','V2560A','V2560A-pre',[msdir2,'/H1.mzXML'],100e-9,10e-6,
       'A2','V2560B','V2560B-pre',[msdir2,'/A2.mzXML'],100e-9,10e-6,
       'E2','V256A-E6','',[msdir2,'/E2.mzXML'],25e-9,10e-6,
       'F2','V256A-F6','',[msdir2,'/F2.mzXML'],25e-9,10e-6,
       'G2','V256A-G6','',[msdir2,'/G2.mzXML'],25e-9,10e-6,
       'H2','V256A-H6','',[msdir2,'/H2.mzXML'],25e-9,10e-6,
       'A1','DMSO','',[msdir3,'/A1.mzXML'],0,20e-6,
       'B1','CDIQ165-B15','',[msdir3,'/B1.mzXML'],100e-9,20e-6,
       'C1','CDIQ245-J11','',[msdir3,'/C1.mzXML'],100e-9,20e-6,
       'D1','CDIQ245-F13','',[msdir3,'/D1.mzXML'],100e-9,20e-6,
       'E1','CDIQ245-P07','',[msdir3,'/E1.mzXML'],100e-9,20e-6,
       'F1','CDIQ485-N16','',[msdir3,'/F1.mzXML'],100e-9,20e-6,
       'G1','V256A-A11','',[msdir3,'/G1.mzXML'],50e-9,20e-6,
       'H1','V5120','',[msdir3,'/H1.mzXML'],100e-9,20e-6,
       'A2','V256A-A1','',[msdir3,'/A2.mzXML'],50e-9,20e-6,
       'B2','V256A-B1','',[msdir3,'/B2.mzXML'],50e-9,20e-6,
       'C2','V256A-C1','',[msdir3,'/C2.mzXML'],50e-9,20e-6,
       'D2','V256A-D1','',[msdir3,'/D2.mzXML'],50e-9,20e-6,
       'E2','V256A-E1','',[msdir3,'/E2.mzXML'],50e-9,20e-6,
       'F2','V256A-F1','',[msdir3,'/F2.mzXML'],50e-9,20e-6,
       'G2','V256A-G1','',[msdir3,'/G2.mzXML'],50e-9,20e-6,
       'H2','V256A-H1','',[msdir3,'/H2.mzXML'],50e-9,20e-6,
       'A3','V256A-A2','',[msdir3,'/A3.mzXML'],50e-9,20e-6,
       'B3','V256A-B2','',[msdir3,'/B3.mzXML'],50e-9,20e-6,
       'C3','V256A-C2','',[msdir3,'/C3.mzXML'],50e-9,20e-6,
       'D3','V256A-D2','',[msdir3,'/D3.mzXML'],50e-9,20e-6,
       'E3','V256A-E2','',[msdir3,'/E3.mzXML'],50e-9,20e-6,
       'F3','V256A-F2','',[msdir3,'/F3.mzXML'],50e-9,20e-6,
       'G3','V256A-G2','',[msdir3,'/G3.mzXML'],50e-9,20e-6,
       'H3','V256A-H2','',[msdir3,'/H3.mzXML'],50e-9,20e-6,
       'A4','V256A-A7','',[msdir3,'/A4.mzXML'],50e-9,20e-6,
       'B4','V256A-B7','',[msdir3,'/B4.mzXML'],50e-9,20e-6,
       'C4','V256A-C7','',[msdir3,'/C4.mzXML'],50e-9,20e-6,
       'D4','V256A-D7','',[msdir3,'/D4.mzXML'],50e-9,20e-6,
       'E4','V256A-E7','',[msdir3,'/E4.mzXML'],50e-9,20e-6,
       'F4','V256A-F7','',[msdir3,'/F4.mzXML'],50e-9,20e-6,
       'G4','V256A-G7','',[msdir3,'/G4.mzXML'],50e-9,20e-6,
       'H4','V256A-H7','',[msdir3,'/H4.mzXML'],50e-9,20e-6,
       'A5','V256A-A8','',[msdir3,'/A5.mzXML'],50e-9,20e-6,
       'B5','V256A-B8','',[msdir3,'/B5.mzXML'],50e-9,20e-6,
       'C5','V256A-C8','',[msdir3,'/C5.mzXML'],50e-9,20e-6,
       'D5','V256A-D8','',[msdir3,'/D5.mzXML'],50e-9,20e-6,
       'E5','CDIQ005-L12','',[msdir3,'/E5.mzXML'],100e-9,20e-6,
       'F5','CDIQ005-N22','',[msdir3,'/F5.mzXML'],100e-9,20e-6,
       'G5','CDIQ045-M11','',[msdir3,'/G5.mzXML'],100e-9,20e-6,
       'H5','CDIQ045-G14','',[msdir3,'/H5.mzXML'],100e-9,20e-6,
       'A6','CDIQ045-N18','',[msdir3,'/A6.mzXML'],100e-9,20e-6,
       'B6','CDIQ125-C17','',[msdir3,'/B6.mzXML'],100e-9,20e-6,
       'C6','CDIQ125-E17','',[msdir3,'/C6.mzXML'],100e-9,20e-6,
       'D6','CDIQ125-J17','',[msdir3,'/D6.mzXML'],100e-9,20e-6,
       'E6','CDIQ125-K21','',[msdir3,'/E6.mzXML'],100e-9,20e-6,
       'F6','CDIQ125-N16','',[msdir3,'/F6.mzXML'],100e-9,20e-6,
       'G6','CDIQ125-O21','',[msdir3,'/G6.mzXML'],100e-9,20e-6,
       'H6','CDIQ165-N09','',[msdir3,'/H6.mzXML'],100e-9,20e-6,
       'A7','CDIQ165-P09','',[msdir3,'/A7.mzXML'],100e-9,20e-6,
       'B7','CDIQ205-G18','',[msdir3,'/B7.mzXML'],100e-9,20e-6,
       'C7','CDIQ205-I19','',[msdir3,'/C7.mzXML'],100e-9,20e-6,
       'D7','CDIQ245-H03','',[msdir3,'/D7.mzXML'],100e-9,20e-6,
       'E7','CDIQ325-C20','',[msdir3,'/E7.mzXML'],100e-9,20e-6,
       'F7','CDIQ325-O09','',[msdir3,'/F7.mzXML'],100e-9,20e-6,
       'G7','CDIQ405-G17','',[msdir3,'/G7.mzXML'],100e-9,20e-6,
       'H7','CDIQ405-M17','',[msdir3,'/H7.mzXML'],100e-9,20e-6,
       'A8','DMSO','',[msdir3,'/A8.mzXML'],0,20e-6
          };
% msdir2 V256A vectors happened to be on same well numbers in QTOF plate
wells={'A3','B3','C3','D3','E3','F3','G3','H3',...
       'A4','B4','C4','D4','E4','F4','G4','H4',...
       'A5','B5','C5','D5','E5','F5','G5','H5',...
       'A6','B6','C6','D6'};
for i=1:length(wells)
  data(end+1,:)={wells{i},['V256A-',wells{i}],'',[msdir2,'/',wells{i},'.mzXML'],25e-9,10e-6};
end
data

data=cell2struct(data,{'well','mixture','name','filename','conc','inject'},2);
for i=1:size(data,1)
  data(i).mzmap=[113,113;500,500];
  data(i).timemap=[0,0;1,1];
  if isempty(data(i).name)
    data(i).name=data(i).mixture;
  end
end

allmz=[];names={};
for i=1:length(qsetup.ADDUCTS)
  for j=1:length(qsetup.mass)
    allmz(end+1)=qsetup.mztarget(j,i);
    names{end+1}=[qsetup.names{j},'[',qsetup.ADDUCTS(i).name,']'];
  end
end

% Load data
if ~exist('mzdata','var')
  mzdata={};
end

for i=1:size(data,1)
  path=data(i).filename;
  if i>length(mzdata) || isempty(mzdata{i}) || ~strcmp(mzdata{i}.path,path)
    % See if we have a .mat file first
    matpath=strrep(path,'.mzXML','.mat');
    if exist(matpath,'file')
      fprintf('Load matfile from %s\n', matpath);
      mztmp=load(matpath);
      fn=fieldnames(mztmp);
      assert(length(fn)==1);  % Saved as a variable
      mzdata{i}=mztmp.(fn{1});
    else
      fprintf('Load mzXML from %s\n', path);
      mztmp=MassSpec(path);
      mztmp.setLoad(data(i).conc*data(i).inject);
      % Prune out some data
      mztmp.filter([0,mztmp.time(end)-5],[min(allmz)-1,max(allmz)+1]);
                                                      %mztmp.keepmz(allmz,'mztol',mztol*2);
      mzdata{i}=mztmp;
      save(matpath,'mztmp');
      fprintf('Saved matfile in %s\n', matpath);
    end
  end
  mzdata{i}.name=data(i).name;   % Reset name to our list
  if ismember(mzdata{i}.path,qsetup.files)
    fprintf('Skipping reload of %s\n', mzdata{i}.name);
    continue;
  end

  n=strsplit(data(i).name,'-');
  % Save record in qsetup
  qsetup.addMS(mzdata{i},data(i).mixture,'group',n{1},'map',struct('mz',data(i).mzmap,'time',data(i).timemap),'sample',data(i).name);
  assert(strcmp(qsetup.files{i},mzdata{i}.path));
end
