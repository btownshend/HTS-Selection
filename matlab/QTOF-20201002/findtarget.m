% Scan a mass spec file for a particular m/z include adduct and isotope search
function findtarget(obj,mzd,id,varargin)
  defaults=struct('mztol',obj.MZFUZZ,'timetol',obj.TIMEFUZZ,'noise',500,'adducts',obj.ADDUCTS);
  args=processargs(defaults,varargin);

  masslist=obj.mass(id)+[obj.ADDUCTS.mass];
  namelist={obj.ADDUCTS.name};
  nadducts=length(masslist);

  isotopes=getisotopes(obj.sdf.sdf(id).Formula,'minabundance',.001);
  tstep=median(diff(mzd.time));
  for i=1:length(args.adducts)
    a=args.adducts(i);
    mz=isotopes(1).mass+a.mass;
    fprintf('%d %s[%s] m/z=%.4f\n', id, obj.names{id}, a.name,obj.mass(id)+a.mass);
    fl=mzd.targetedFeatureDetect(mz,'names',{a.name},'mztol',args.mztol,'noise',args.noise,'timetol',args.timetol);
    fld=fl.deconvolve('noise',500);
    for j=1:length(fld.features)
      f=fld.features(j);
      fprintf(' %-41.41s %.4f d=%3.0f T=%5.2f, I=%7.0f\n', sprintf('%s@%.2f',isotopes(1).name,f.time), f.mz, (f.mz-mz)*1e4, f.time, f.intensity);
      %fprintf('  timerange=[%.3f,%.3f]\n', f.timerange);
      for k=2:length(isotopes)
        if isotopes(k).abundance*f.intensity < args.noise/2
          continue;
        end
        mzi=isotopes(k).mass+a.mass;
        fli=mzd.targetedFeatureDetect(mzi,'names',{isotopes(k).name},'mztol',args.mztol,'noise',0,'timetol',tstep/2,'rt',f.time);
        fprintf('  %-40.40s %.4f',isotopes(k).name,mzi);
        if length(fli.features)>0
          fi=fli.features(1);
          assert(fi.time==f.time);
          assert(fi.npeaks==1);
          fprintf(' d=%3.0f T=%5.2f, I=%7.0f Rel=%4.1f%%',(fi.mz-mzi)*1e4, fi.time, fi.intensity,fi.intensity/f.intensity*100);
        else
          fprintf(' d=%3.0f T=%5.2f, I=%7.0f Rel=%4.1f%%',0, f.time, 0,0);
        end
        fprintf(' Expected %4.1f%%',isotopes(k).abundance/isotopes(1).abundance*100);
        fprintf('\n');
      end
    end
  end
end



