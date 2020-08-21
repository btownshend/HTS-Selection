trace=300;
compounds.assignTimes('clear',true,'timetol',0.35,'mztol',40e-4,'trace',trace);
compounds.checksensitivity(1);
compounds.assignTimes('clear',false,'timetol',0.35,'mztol',40e-4,'trace',trace);
compounds.assignTimes('clear',false,'timetol',0.35,'normicrange',[0.5,4.5],'trace',trace);
compounds.assignTimes('clear',false,'timetol',0.35,'maxFP',1,'trace',trace);
compounds.assignTimes('clear',false,'timetol',0.35,'maxFN',1,'trace',trace);
compounds.assignTimes('clear',false,'timetol',0.35,'maxFN',1,'maxFP',1,'trace',trace);
compounds.checksensitivity(1);
compounds.assignTimes('clear',false,'timetol',0.35,'maxFN',2,'maxFP',2,'trace',trace);

