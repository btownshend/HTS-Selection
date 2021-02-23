-- MySQL Script generated by MySQL Workbench
-- Thu Feb 18 22:35:08 2021
-- Model: New Model    Version: 1.0
-- MySQL Workbench Forward Engineering

SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='TRADITIONAL,ALLOW_INVALID_DATES';

-- -----------------------------------------------------
-- Schema compounds
-- -----------------------------------------------------

-- -----------------------------------------------------
-- Schema compounds
-- -----------------------------------------------------
CREATE SCHEMA IF NOT EXISTS `compounds` DEFAULT CHARACTER SET utf8 ;
USE `compounds` ;

-- -----------------------------------------------------
-- Table `compounds`.`compounds`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `compounds`.`compounds` (
  `compound` INT NOT NULL,
  `inChiKey` CHAR(27) NULL,
  `monoisotopicMass` FLOAT NULL,
  `formula` VARCHAR(200) NULL,
  PRIMARY KEY (`compound`),
  UNIQUE INDEX `inChiKey_UNIQUE` (`inChiKey` ASC))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `compounds`.`mixtures`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `compounds`.`mixtures` (
  `mixture` INT NOT NULL,
  `name` VARCHAR(45) NULL,
  `solvent` INT NOT NULL,
  PRIMARY KEY (`mixture`),
  INDEX `fk_Samples_Compounds1_idx` (`solvent` ASC),
  CONSTRAINT `fk_Samples_Compounds1`
    FOREIGN KEY (`solvent`)
    REFERENCES `compounds`.`compounds` (`compound`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `compounds`.`contents`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `compounds`.`contents` (
  `mixture` INT NOT NULL,
  `compound` INT NOT NULL,
  `concentration` FLOAT NOT NULL,
  PRIMARY KEY (`mixture`, `compound`),
  INDEX `fk_Samples_has_Compounds_Compounds1_idx` (`compound` ASC),
  INDEX `fk_Samples_has_Compounds_Samples_idx` (`mixture` ASC),
  CONSTRAINT `fk_Samples_has_Compounds_Samples`
    FOREIGN KEY (`mixture`)
    REFERENCES `compounds`.`mixtures` (`mixture`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_Samples_has_Compounds_Compounds1`
    FOREIGN KEY (`compound`)
    REFERENCES `compounds`.`compounds` (`compound`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `compounds`.`domains`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `compounds`.`domains` (
  `domain` INT NOT NULL,
  `name` VARCHAR(45) NULL,
  PRIMARY KEY (`domain`))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `compounds`.`names`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `compounds`.`names` (
  `compound` INT NOT NULL,
  `domain` INT NOT NULL,
  `name` VARCHAR(45) NOT NULL,
  INDEX `fk_Identifiers_Compounds1_idx` (`compound` ASC),
  PRIMARY KEY (`compound`, `domain`),
  INDEX `fk_Identifiers_idDomains1_idx` (`domain` ASC),
  UNIQUE INDEX `uniqueid` (`domain` ASC, `name` ASC),
  CONSTRAINT `fk_Identifiers_Compounds1`
    FOREIGN KEY (`compound`)
    REFERENCES `compounds`.`compounds` (`compound`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_Identifiers_idDomains1`
    FOREIGN KEY (`domain`)
    REFERENCES `compounds`.`domains` (`domain`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `compounds`.`proptypes`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `compounds`.`proptypes` (
  `proptype` INT NOT NULL,
  `name` VARCHAR(45) NOT NULL,
  PRIMARY KEY (`proptype`))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `compounds`.`properties`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `compounds`.`properties` (
  `compound` INT NOT NULL,
  `proptype` INT NOT NULL,
  `value` VARCHAR(45) NOT NULL,
  INDEX `fk_compoundProps_Compounds1_idx` (`compound` ASC),
  INDEX `fk_compoundProps_properties1_idx` (`proptype` ASC),
  PRIMARY KEY (`proptype`, `compound`),
  CONSTRAINT `fk_compoundProps_Compounds1`
    FOREIGN KEY (`compound`)
    REFERENCES `compounds`.`compounds` (`compound`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_compoundProps_properties1`
    FOREIGN KEY (`proptype`)
    REFERENCES `compounds`.`proptypes` (`proptype`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


SET SQL_MODE=@OLD_SQL_MODE;
SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS;
SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;