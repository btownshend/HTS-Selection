CREATE DATABASE hts;

USE hts;

DROP TABLE mscontains;
DROP TABLE isopeaks;
DROP TABLE features;
DROP TABLE msruns;
DROP TABLE isotopes;
DROP TABLE compounds;
DROP TABLE formulas;
DROP TABLE adducts;


CREATE TABLE formulas (
       formula_pk INTEGER AUTO_INCREMENT,  PRIMARY KEY(formula_pk),
       formula VARCHAR(100) NOT NULL
);

CREATE TABLE isotopes (
       formula_pk INTEGER NOT NULL,  FOREIGN KEY(formula_pk) REFERENCES formulas(formula_pk),
       isotope INTEGER AUTO_INCREMENT,  PRIMARY KEY(isotope),
       formula VARCHAR(100) NOT NULL,
       mass FLOAT NOT NULL,  -- exact mass
       abundance FLOAT NOT NULL  -- abundance of this isotope (not relative to monoisotopic mass)
);

CREATE TABLE compounds (
       compound INTEGER AUTO_INCREMENT,  PRIMARY KEY(compound),
       id VARCHAR(20) NOT NULL,   -- e.g. 'CDIV-604E11'
       formula_pk INTEGER NOT NULL, FOREIGN KEY(formula_pk) REFERENCES formulas(formula_pk),
       unique(id)
);

CREATE TABLE msruns (
       msrun INTEGER AUTO_INCREMENT, PRIMARY KEY(msrun),
       created DATE,
       name VARCHAR(50) NOT NULL,
       filename VARCHAR(50) NOT NULL
);	

CREATE TABLE mscontains (
       msrun INTEGER, FOREIGN KEY(msrun) REFERENCES msruns(msrun),
       compound INTEGER, FOREIGN KEY(compound) REFERENCES compounds(compound),
       PRIMARY KEY(msrun,compound)
);

CREATE TABLE adducts (
       adduct INTEGER AUTO_INCREMENT, PRIMARY KEY(adduct),
       name VARCHAR(20) NOT NULL,
       mass FLOAT NOT NULL
);

-- a feature is the isotope pattern of a particular compound(formula) at one elution time
CREATE TABLE features (
	   feature INTEGER AUTO_INCREMENT, PRIMARY KEY(feature),
       msrun INTEGER, FOREIGN KEY(msrun) REFERENCES msruns(msrun),
       formula_pk INTEGER NOT NULL, FOREIGN KEY(formula_pk) REFERENCES formulas(formula_pk),
       adduct INTEGER, FOREIGN KEY(adduct) REFERENCES adducts(adduct),
       rt FLOAT NOT NULL,  -- exact RT of peak interrogated for isotopes
       mztol FLOAT,  -- mztol used to find isotopes
       UNIQUE(msrun,formula_pk,adduct,rt)
);
    
CREATE TABLE isopeaks (
       feature INTEGER, FOREIGN KEY(feature) REFERENCES features(feature) ON DELETE CASCADE,
       isotope INTEGER, FOREIGN KEY(isotope) REFERENCES isotopes(isotope),
       obsmz FLOAT,  -- observed m/z if peak found
       ioncount FLOAT NOT NULL, 
       PRIMARY KEY(feature, isotope)
);
   
create or replace view v_features as
SELECT feature,msrun,f.formula,a.name,a.mass amass,rt,mztol
FROM features fe, formulas f, adducts a
WHERE fe.formula_pk=f.formula_pk
AND a.adduct=fe.adduct;

create or replace view v_isopeaks as
select fe.msrun,ip.feature,fe.rt,i.formula,fe.amass+i.mass expmz,obsmz,i.abundance,ioncount
FROM isopeaks ip, isotopes i, v_features fe
WHERE ip.isotope=i.isotope
AND fe.feature=ip.feature;
