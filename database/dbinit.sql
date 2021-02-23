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

-- Make all the CDIQ mixtures use 4-digit plate number
UPDATE mixtures
SET name=REPLACE(name,'CDIQ','CDIQ0')
WHERE name LIKE 'CDIQ%'
AND LENGTH(name)==11;

-- Add CDIV0565-{row} mixtures
INSERT INTO mixtures(name,solvent) 
SELECT CONCAT('CDIV0565-',LEFT(RIGHT(n.name,3),1)),sv.compound	
FROM compounds sv, names n
WHERE sv.name='DMSO'
AND n.domain=(SELECT domain FROM domains d WHERE d.name='CDIV') 
AND n.name LIKE '565%02';

INSERT INTO contents(mixture,compound,concentration) 
SELECT m.mixture,n.compound, 5e-3/10
FROM mixtures m, names n
WHERE n.domain=(SELECT domain FROM domains d WHERE d.name='CDIV') 
AND n.name LIKE CONCAT('565',RIGHT(m.name,1),'%')
AND m.name LIKE 'CDIV0565-_';

-- Empty solution (DMSO or Water)
INSERT INTO mixtures(name,solvent) 
SELECT sv.name,sv.compound	
FROM compounds sv
WHERE sv.name IN ('Water','DMSO');


-- TODO:
-- Should add entries in mixtures for CDIV960 (S6), including column, row, plate, etc sums
      	  -- see NGS88 and prior
-- Should add entries in mixtures for Theo, etc
