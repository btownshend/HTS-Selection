compounds.assignTimes('clear',true,'timetol',0.35,'mztol',40e-4);
compounds.checksensitivity(1);
compounds.assignTimes('clear',false,'timetol',0.35,'mztol',40e-4);
compounds.assignTimes('clear',false,'timetol',0.35,'normicrange',[0.5,4.5]);
compounds.assignTimes('clear',false,'timetol',0.35,'maxFP',1);
compounds.assignTimes('clear',false,'timetol',0.35,'maxFN',1);
compounds.assignTimes('clear',false,'timetol',0.35,'maxFN',1,'maxFP',1);
compounds.checksensitivity(1);

