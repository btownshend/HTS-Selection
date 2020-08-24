% Process contaminants database
datadir='../../data/MassSpec/MaConDa';
infile=[datadir,'/MaConDa__v1_0__extensive.csv'];
outfile=[datadir,'/contaminants.csv'];
x=readtable(infile);
sel=strcmp(x.ion_mode,'POS') & x.exact_adduct_mass~=0;
% Only use ID, name,formula,m/z,RT
x.RT(:)=0;
for i=1:size(x,1)
  x.fullname{i}=[x.name{i},x.ion_form{i}];
end
keepvars={'id','fullname','formula','exact_adduct_mass','RT'};
indices=cellfun(@(z) find(strcmp(x.Properties.VariableNames,z)),keepvars);
xout=x(sel,indices);
writetable(xout,outfile);

