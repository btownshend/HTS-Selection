fsel=[];
for col=1:8
  for row='A':'H'
    if col==8 && row>='E'
      break;
    end
    vname=sprintf('V256A-%c%d',row,col);
    ind=find(strcmp(qsetup.samples,vname));
    if isempty(ind)
      fprintf('Missing vector: %s\n', vname);
    else
      fsel(end+1)=ind;
    end
  end
end

    
v256=qsetup.copypart('fsel',fsel);
timetol=0.2;
mztol=.0005;
exphits=3;
trace=[];

v256.fsens(:)=1;
v256.assignTimes('clear',true,'timetol',timetol,'mztol',mztol,'normicrange',[0,1000],'trace',trace,'minhits',exphits);
v256.checksensitivity();  % Initial fsens
v256.assignTimes('clear',true,'timetol',timetol,'mztol',mztol,'normicrange',[0.05,20],'trace',trace,'minhits',exphits);
v256.checksensitivity();
v256.assignTimes('clear',false,'timetol',timetol,'mztol',mztol,'trace',trace,'minhits',exphits);
v256.assignTimes('clear',false,'timetol',timetol,'mztol',mztol,'trace',trace,'minhits',exphits-1);
v256.checksensitivity();
for fn=1:3
  v256.assignTimes('clear',false,'timetol',timetol,'maxFP',fn,'mztol',mztol,'trace',trace,'minhits',exphits-1);
end
v256.checksensitivity();
v256.platesummary;
save('v256.mat','v256');

report=v256.report();
ds=datestr(now,'YYYYmmDDHHMM');
resultsdir='../../results/';
writetable(report,'/tmp/report.csv','QuoteStrings',true);
system(sprintf('sed -e "s/NaN//g" /tmp/report.csv > %s',[resultsdir,sprintf('report-%s.csv',ds)]));
