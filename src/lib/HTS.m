classdef HTS < handle
  properties
    molecules;
  end
  
  methods
    function obj=HTS()
      obj.molecules=[];
    end
    
    function load(obj, file)
    % Load data from Excel file from HTBC
      fd=fopen(file,'r');
      if fd<0
        error('Unable to open %s', file);
      end
      % IDNUMBER,IDNUMBER,compound corp reg number,BatchMW,Plate,Well,SALTDATA,mol weight *structure,logP,logD,logSw,H acceptor,H donor,B rotN,N O,Amount,
      data=textscan(fd,'%s %s %s %f %s %s %s %f %f %f %f %f %f %f %f %f','HeaderLines',1,'Delimiter',',');
      fprintf('Loaded %d entries from %s\n', length(data{1}), file);
      d=struct('id1',data{1},'id2',data{2},'ccrnum',data{3},'batchMW',num2cell(data{4}),'plate',data{5},'well',data{6},'saltdata',data{7},'structMW',num2cell(data{8}),'logP',num2cell(data{9}),'logD',num2cell(data{10}),'logSw',num2cell(data{11}),'Hacceptor',num2cell(data{12}),'Hdonor',num2cell(data{13}),'BrotN',num2cell(data{14}),'NO',num2cell(data{15}),'Amount',num2cell(data{16}));
      obj.molecules=[obj.molecules;d];
      fclose(fd);
    end

    function mw=getmw(obj,plate)
    % Get MW's from given plate
      mw=nan(8,12);
      onplate=strcmp({obj.molecules.plate},plate);
      if sum(onplate)==0
        error('No entries for plate %s', plate);
      end
      %fprintf('Have %d entries for plate %s\n', plate);
      %{obj.molecules(onplate).well}
      mp=obj.molecules(onplate);
      for r=1:8
        for c=2:11
          well=sprintf('%c%02d',char(r+'A'-1),c);
          sel=strcmp({mp.well},well);
          if sum(sel)==1
            mw(r,c)=mp(sel).structMW;
          elseif sum(sel)~=0
            error('Have %d entries for plate %s, well %s', plate, well);
          end
        end
      end
      
    end
  end % methods
end % classdef



