% Load data from Excel file from HTBC
file='SET_50000.csv';
fd=fopen(file,'r');
if fd<0
  error('Unable to open %s', file);
end
% IDNUMBER,IDNUMBER,compound corp reg number,BatchMW,Plate,Well,SALTDATA,mol weight *structure,logP,logD,logSw,H acceptor,H donor,B rotN,N O,Amount,
data=textscan(fd,'%s %s %s %f %s %s %s %f %f %f %f %f %f %f %f %f','HeaderLines',1,'Delimiter',',');
d=struct('id1',data{1},'id2',data{2},'ccrnum',data{3},'batchMW',num2cell(data{4}),'plate',data{5},'well',data{6},'saltdata',data{7},'structMW',num2cell(data{8}),'logP',num2cell(data{9}),'logD',num2cell(data{10}),'logSw',num2cell(data{11}),'Hacceptor',num2cell(data{12}),'Hdonor',num2cell(data{13}),'BrotN',num2cell(data{14}),'NO',num2cell(data{15}),'Amount',num2cell(data{16}));
fclose(fd);



