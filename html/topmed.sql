/* --------------------------------------------------------------
#       NHLBI TopMed tracking entries
#--------------------------------------------------------------- 
create database nhlbi;
GRANT ALL    ON nhlbi.* TO sqlnhlbi@"localhost"   IDENTIFIED BY 'password';
GRANT SELECT ON nhlbi.* TO sqlnhlbiro@"localhost" IDENTIFIED BY 'password';
flush privileges;

#   When you see instruction to do:  mysqladmin flush-hosts do:
ssh d
ssh root@g
mysql -h f-db -p
flush hosts;
#   Verify it works by connecting on host that is having the problem

#   Clone the old database to the new
mysqldump -u sqlnhlbi --password=pw -h f-db nhlbi > /tmp/nhlbi.sql
mysql -u sqlnhlbi --password=pw -h localhost nhlbi < /tmp/nhlbi.sql
*/

/* This table is used to control when/if an operation is permitted
   If an entry is in this table, it means it is disabled */
DROP TABLE IF EXISTS permissions;
CREATE TABLE permissions (
  id           INT         NOT NULL AUTO_INCREMENT,
  datayear     INT         NOT NULL,            # 0 = all
  centerid     INT         NOT NULL,            # 0 = all
  centername   VARCHAR(16) NOT NULL,
  runid        INT         NOT NULL,            # 0 = all
  dirname      VARCHAR(64) NOT NULL,
  operation    CHAR(12),
  PRIMARY KEY  (id)
);

/* Lists all the centers we get BAMs from.
    Datamethod is push or pull depending on how data gets to us. */
DROP TABLE IF EXISTS centers;
CREATE TABLE centers (
  centerid     INT NOT NULL AUTO_INCREMENT, /* Referenced in nhlbi_qc_metrics */
  centername   VARCHAR(16) NOT NULL,        /* Referenced in nhlbi_qc_metrics */
  datamethod   VARCHAR(8) NOT NULL DEFAULT 'pull',
  centerdesc   VARCHAR(96) NOT NULL,
  designdesc   TEXT NOT NULL,
  PRIMARY KEY  (centerid)
);
INSERT INTO centers (centername,centerdesc,designdesc) VALUES('broad', 'Broad Institute','Illumina sequencing of Homo sapiens via random selection');
INSERT INTO centers (centername,centerdesc,designdesc) VALUES('illumina', '30x whole genome sequencing using Illumina TruSeq PCR-free library protocol with 500ng-1ug input DNA');
INSERT INTO centers (centername,centerdesc,designdesc) VALUES('nygc', 'New York Genome Center','Whole genome sequencing using Illumina TruSeq PCR-free DNA library preparation with 500ng input DNA, sequenced to >30x mean coverage with 2x150bp reads on HiSeq X.');
INSERT INTO centers (centername,centerdesc,designdesc) VALUES('uw', '30x Illumina whole genome sequencing using Kapa PCR-free DNA library preparation with 500ng input DNA');
INSERT INTO centers (centername,centerdesc,designdesc) VALUES('baylor', 'BCM','30x Illumina whole genome sequencing using Kapa PCR-free DNA library preparation with 500ng input DNA');
INSERT INTO centers (centername,centerdesc,designdesc) VALUES('macrogen', 'Macrogen','30x Illumina whole genome sequencing using Kapa PCR-free DNA library preparation with 500ng input DNA');
INSERT INTO centers (centername,centerdesc,designdesc) VALUES('washu', 'McDonnell Genome Center at Washington University','30x Illumina whole genome sequencing using Kapa PCR-free DNA library preparation with 500ng input DNA');

/*  The second of these is used to overwhelm a bug on the part of NCBI */
UPDATE centers set centerdesc='New York Genome Center' where centerid=3;
UPDATE centers set centerdesc='NYGC' where centerid=3;

/* Lists all runs (directories of data). Must be tied to a center */
DROP TABLE IF EXISTS runs;
CREATE TABLE runs (
  runid        INT NOT NULL AUTO_INCREMENT,     /* Referenced in nhlbi_qc_metrics */
  centerid     INT,
  dirname      VARCHAR(64) NOT NULL,
  status       VARCHAR(256),
  count        INT,
  datayear     INT DEFAULT 3,          /* Year of project: 1, 2 ... */
  offsite      CHAR(1) DEFAULT 'N',    /* Original files kept offsite (N,Y,D) */
  xmlfound     INT DEFAULT 0,          /* Remove this */
  arrived      CHAR(1) DEFAULT 'N',    /* Y or N that all files arrived for this run */

  dateinit     VARCHAR(12),
  datecomplete VARCHAR(12),

  comments     TEXT,
  PRIMARY KEY  (runid)
);
CREATE INDEX index_centerid_dirname ON runs(centerid,dirname);
CREATE INDEX index_dirname ON runs(dirname);

 /* All BAMs we know about (e.g. NWDID)  Must be tied to a run */
DROP TABLE IF EXISTS bamfiles;
CREATE TABLE bamfiles (
  bamid        INT         NOT NULL AUTO_INCREMENT, /* Referenced in nhlbi_qc_metrics */
  runid        INT         NOT NULL,
  bamname      VARCHAR(96) NOT NULL,
  base_coord   INT DEFAULT 0,
  library_name VARCHAR(96),
  nominal_sdev INT DEFAULT 0,
  nominal_length  INT DEFAULT 0,
  cramname     VARCHAR(96) NOT NULL,
  cramchecksum VARCHAR(96) NOT NULL,
  b37cramchecksum VARCHAR(96) NOT NULL,
  b38cramchecksum VARCHAR(96) NOT NULL,
  b38craichecksum VARCHAR(96) NOT NULL,
  bamflagstat  BIGINT UNSIGNED DEFAULT NULL,
  cramflagstat BIGINT UNSIGNED DEFAULT NULL,
  b37flagstat  BIGINT UNSIGNED DEFAULT NULL,
  b38flagstat  BIGINT UNSIGNED DEFAULT NULL,
  datemapping_b38  datetime DEFAULT NULL,   /* Referenced in nhlbi_qc_metrics */
  datemapping_b37  datetime DEFAULT NULL,   /* Referenced in nhlbi_qc_metrics */
  studyname    VARCHAR(96) NOT NULL,    /* Referenced in nhlbi_qc_metrics */
  piname       VARCHAR(96),             /* Referenced in nhlbi_qc_metrics */
  phs          VARCHAR(12),
  phs_consent_short_name VARCHAR(24),
  phs_sra_sample_id VARCHAR(24),
  phs_sra_data_details VARCHAR(255),
  emsg         VARCHAR(255),
  checksum     VARCHAR(96) NOT NULL,
  expt_sampleid VARCHAR(24),            /* Referenced in nhlbi_qc_dometrics, NWDID */
  nwdid_known  CHAR(1) DEFAULT 'N',     /* Sample is known to NCBI */
  poorquality  CHAR(1) DEFAULT 'N',     /* Quality of sample is poor, do not use */
  send2aws     CHAR(1) DEFAULT 'Y',     /* Send this sample to AWS */

  datearrived  VARCHAR(12),             /* Referenced in nhlbi_qc_metrics */
  bamsize      VARCHAR(16) DEFAULT 0,   /* Referenced in nhlbi_qc_metrics */
  datayear     INT DEFAULT 3,           /* Year of project: 1, 2 ... */
  build        VARCHAR(4) DEFAULT '38', /* Build original input file user, 37, 38 etc */
  dateinit     VARCHAR(12),             /* Referenced in nhlbi_qc_metrics */
  bamname_orig VARCHAR(96) NOT NULL,
  offsite      CHAR(1) DEFAULT 'N',     /* Original files kept offsite (N,Y,D) */

/* Fields to track state for each step */
/*
my $NOTSET    = 0;            # Not set
my $REQUESTED = 1;            # Task requested
my $SUBMITTED = 2;            # Task submitted to be run
my $STARTED   = 3;            # Task started
my $DELIVERED = 19;           # Data delivered, but not confirmed
my $COMPLETED = 20;           # Task completed successfully
my $IGNORETHIS = 80;          # Task is to be ignored
my $FAILEDCHECKSUM = 88;      # Task failed, because checksum at NCBI bad
my $CANCELLED = 89;           # Task cancelled
my $FAILED    = 99;           # Task failed
*/
  state_arrive   INT DEFAULT 0,
  state_verify   INT DEFAULT 0,
  state_gcebackup INT DEFAULT 0,
  state_cram     INT DEFAULT 0,
  state_qplot    INT DEFAULT 0,     /* Referenced in nhlbi_qc_metrics */
  state_b37      INT DEFAULT 0,     /* Referenced in nhlbi_qc_metrics */

  state_gce38push  INT DEFAULT 0,   /* Handling b38 data in Google Cloud */
  state_gce38pull  INT DEFAULT 0,
  state_b38      INT DEFAULT 0,     /* Referenced in nhlbi_qc_metrics */

  state_gce38bcf  INT DEFAULT 0,    /* Created BCF data */
  state_gce38cpbcf  INT DEFAULT 0,  /* Copy BCF data to GCE BCF bucket */

  state_gce38copy INT DEFAULT 0,    /* Copy local data to SHARE bucket */
  state_aws38copy INT DEFAULT 0,    /* Copy local data to AWS bucket */

  state_ncbiexpt INT DEFAULT 0,     /* Year one, experiment defined at NCBI (X) */
  state_ncbiorig INT DEFAULT 0,     /* Original input file as bam or cram (S) */
  state_ncbib37  INT DEFAULT 0,     /* Remapped cram build 37 as cram (P) */

  state_fix INT DEFAULT 0,          /* Track efforts to fix screwups */

  PRIMARY KEY  (bamid)
);

/*   Remove this from bamfiles
     Remove this from runs
*/


/*   Handy queries
    ALTER TABLE bamfiles ADD COLUMN datebai VARCHAR(12) AFTER datebackup;
    ALTER TABLE runs ADD  COLUMN offsitebackup CHAR(1) DEFAULT 'N' AFTER datayear;
    ALTER TABLE bamfiles CHANGE expt_sampleid expt_sampleid VARCHAR(24);
    
    ALTER TABLE bamfiles DROP COLUMN colname;
*/
 
CREATE INDEX index_runid   ON bamfiles(runid);
CREATE INDEX index_nwdid   ON bamfiles(expt_sampleid);
CREATE INDEX index_refname ON bamfiles(refname);
CREATE INDEX index_datayear ON bamfiles(datayear);
CREATE INDEX index_bamflagstat ON bamfiles(bamflagstat);
CREATE INDEX index_cramflagstat ON bamfiles(cramflagstat);
CREATE INDEX index_b37flagstat ON bamfiles(b37flagstat);
CREATE INDEX index_b38flagstat ON bamfiles(b38flagstat);
CREATE INDEX index_studyname ON bamfiles(studyname);
CREATE INDEX index_piname ON bamfiles(piname);
CREATE INDEX index_state_verify ON bamfiles(state_verify);
CREATE INDEX index_state_cram ON bamfiles(state_cram);
CREATE INDEX index_state_bai ON bamfiles(state_bai);
CREATE INDEX index_state_qplot ON bamfiles(state_qplot);
CREATE INDEX index_state_b37 ON bamfiles(state_b37);
CREATE INDEX index_state_b38 ON bamfiles(state_b38);
CREATE INDEX index_state_gce38push ON bamfiles(state_gce38push);
CREATE INDEX index_state_gce38pull ON bamfiles(state_gce38pull);
CREATE INDEX index_state_gce38post ON bamfiles(state_gce38post);
CREATE INDEX index_state_bcf ON bamfiles(state_bcf);
ALTER TABLE bamfiles ADD UNIQUE (expt_sampleid);

-- MySQL Workbench Forward Engineering

SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='TRADITIONAL,ALLOW_INVALID_DATES';

/*  This table is regularly loaded from a tab delimited summary file from NCBI
    The NCBI data is loaded into this table to make it easier for us to garner
    the state of data sent there.

    mysql -h f-db.sph.umich.edu -u USER --password=PASSWORD --local-infile=1 nhlbi    
    load data local infile 'ncbi_summary' into table ncbi_summary;

*/
DROP TABLE IF EXISTS ncbi_summary;
CREATE TABLE ncbi_summary (
  realm         VARCHAR(12),
  upload_id     VARCHAR(12),
  upload_date   VARCHAR(20),
  file_name     VARCHAR(64),
  file_size     VARCHAR(12),
  file_md5sum   VARCHAR(32),
  upload_name   VARCHAR(64),
  upload_size   VARCHAR(24),
  upload_md5sum VARCHAR(32),
  file_status   VARCHAR(20),
  file_type	    VARCHAR(20),
  load_date     VARCHAR(24),
  file_error    VARCHAR(255),
  submissions   VARCHAR(64),
  loaded_runs   VARCHAR(64),
  unloaded_runs VARCHAR(64),
  suppressed_runs   VARCHAR(64),
  loaded_analyses   VARCHAR(64),
  unloaded_analyses VARCHAR(64),
  suppressed_analyses   VARCHAR(64),
  PRIMARY KEY  (file_name)
);


/* ####################################################
   Daily statisitics for steps
   #################################################### */
DROP TABLE IF EXISTS stepstats;
CREATE TABLE stepstats (
  yyyymmdd CHAR(10) NOT NULL,
  count_verify      INT DEFAULT 0,
  avetime_verify    INT DEFAULT 0,
  count_bai         INT DEFAULT 0,
  avetime_bai       INT DEFAULT 0,
  count_qplot       INT DEFAULT 0,
  avetime_qplot     INT DEFAULT 0,
  count_cram        INT DEFAULT 0,
  avetime_cram      INT DEFAULT 0,
  count_expt        INT DEFAULT 0,
  avetime_expt      INT DEFAULT 0,
  ncbicount_expt    INT DEFAULT 0,
  count_orig        INT DEFAULT 0,
  avetime_orig      INT DEFAULT 0,
  ncbicount_orig    INT DEFAULT 0,
  count_b37         INT DEFAULT 0,
  avetime_b37       INT DEFAULT 0,
  ncbicount_b37     INT DEFAULT 0,
  count_gcepush     INT DEFAULT 0,
  avetime_gcepush   INT DEFAULT 0,
  ncbicount_gcepush INT DEFAULT 0,
  count_gcepull     INT DEFAULT 0,
  avetime_gcepull   INT DEFAULT 0,
  ncbicount_gcepull INT DEFAULT 0,
  count_gcecopy     INT DEFAULT 0,
  avetime_gcecopy   INT DEFAULT 0,
  ncbicount_gcecopy INT DEFAULT 0,
  count_gcecpbcf    INT DEFAULT 0,
  avetime_gcecpbcf  INT DEFAULT 0,
  ncbicount_gcecpbcf INT DEFAULT 0,
  count_gcebackup     INT DEFAULT 0,
  avetime_gcebackup   INT DEFAULT 0,
  ncbicount_gcebackup INT DEFAULT 0,
  count_awscopy     INT DEFAULT 0,
  avetime_awscopy   INT DEFAULT 0,
  ncbicount_awscopy INT DEFAULT 0,

  count_fix         INT DEFAULT 0,
  avetime_fix       INT DEFAULT 0,
  ncbicount_fix     INT DEFAULT 0,

  count_b38         INT DEFAULT 0,
  avetime_b38       INT DEFAULT 0,
  ncbicount_b38     INT DEFAULT 0,
  count_bcf         INT DEFAULT 0,
  avetime_bcf       INT DEFAULT 0,
  ncbicount_bcf     INT DEFAULT 0,

  bamcount           INT DEFAULT 0,      /* Count of all arrived bams */
  errcount           INT DEFAULT 0,      /* Count of all errors for bams */
  errorigcount       INT DEFAULT 0,      /* Count of original bams sent to NCBI in error */
  loadedorigbamcount INT DEFAULT 0,      /* Count of loaded original BAMs at NCBI */
  errckorigcount     INT DEFAULT 0,      /* Count of original bams at NCBI with checsum error */
  errb37count        INT DEFAULT 0,      /* Count of primary bams sent to NCBI in error */
  loadedb37bamcount  INT DEFAULT 0,      /* Count of loaded primary BAMs at NCBI */
  errckb37count      INT DEFAULT 0,      /* Count of primary bams at NCBI with checsum error */
  errb38count        INT DEFAULT 0,      /* Count of tertiary bams sent to NCBI in error */
  loadedb38bamcount  INT DEFAULT 0,      /* Count of loaded tertiary BAMs at NCBI */
  errckb38count      INT DEFAULT 0,      /* Count of tertiary bams at NCBI with checsum error */
  PRIMARY KEY  (yyyymmdd)
);

/* ####################################################
   View of Kevin's web interface columns

DROP TABLE IF EXISTS qc_results;
CREATE TABLE qc_results (
  id      INT PRIMARY KEY AUTO_INCREMENT,
  bam_id  int,
  pct_freemix         float,
  n_reads_m           float,
  pct_mapped          float,
  pct_mq0             float,
  pct_paired          float,
  pct_prop_paired     float,
  mapped_gb           float,
  q20_gb              float,
  pct_q20_base        float,
  mean_depth          float,
  pct_genome_cov      float,
  isize_mode          float,
  isize_median        float,
  pct_dups            float,
  pct_genome_dp5      float,
  pct_genome_dp10     float,
  pct_genome_dp20     float,
  pct_genome_dp30     float,
  vmr_depth           float,
  depth_q10           float,
  depth_q20           float,
  depth_q30           float,
  raw_base_gb         float,
  pct_overlap_reads   float,
  pct_overlap_bases   float,
  isize_iqr           float,
  isize_stdev         float,
  gc_depth_0_1        float,
  gc_depth_1_5        float,
  gc_depth_5_25       float,
  gc_depth_25_75      float,
  gc_depth_75_95      float,
  gc_depth_95_99      float,
  gc_depth_99_100     float,
  library_size_m      float,
  qc_fail             int  DEFAULT=0,
  qc_pass             int  DEFAULT=0,
  qc_flagged          int  DEFAULT=0,
  created_at          datetime,                                
  modified_at         timestamp,        
  PRIMARY KEY  (id)
);


This entire VIEW exists because Kevin was porting some messy Hyun code
to make it faster to render and didn't want to completely rewrite
Hyuns dashboard, although it should be at some point. To maintain the
data structure that Hyun had I created the view from the imported qc
data to match what he had in place so Kevin had an easier go of
things. Presumably the conditions that are being tested could be
redefined in the VIEW to get updated results without updating the data
in the qc results table. That said, unless you are redoing the dashboard,
I'd leave this VIEW alone.

CREATE 
    ALGORITHM = UNDEFINED 
    DEFINER = `sqlnhlbi`@`%` 
    SQL SECURITY DEFINER
VIEW `nhlbi_qc_metrics` AS
    SELECT 
        `b`.`expt_sampleid` AS `sample_id`,
        `b`.`piname` AS `pi_name`,
        `b`.`studyname` AS `study`,
        `c`.`centername` AS `center`,
        '2016-04-15' AS `seq_date`,
        FROM_UNIXTIME(`b`.`dateinit`) AS `bam_date`,
        `q`.`pct_freemix` AS `pct_freemix`,
        `q`.`n_reads_m` AS `n_reads_m`,
        `q`.`pct_mapped` AS `pct_mapped`,
        `q`.`pct_mq0` AS `pct_mq0`,
        `q`.`pct_paired` AS `pct_paired`,
        `q`.`pct_prop_paired` AS `pct_prop_paired`,
        `q`.`mapped_gb` AS `mapped_gb`,
        `q`.`q20_gb` AS `q20_gb`,
        `q`.`pct_q20_base` AS `pct_q20_base`,
        `q`.`mean_depth` AS `mean_depth`,
        `q`.`pct_genome_cov` AS `pct_genome_cov`,
        `q`.`isize_mode` AS `isize_mode`,
        `q`.`isize_median` AS `isize_median`,
        `q`.`pct_dups` AS `pct_dups`,
        `q`.`pct_genome_dp5` AS `pct_genome_dp5`,
        `q`.`pct_genome_dp10` AS `pct_genome_dp10`,
        `q`.`pct_genome_dp20` AS `pct_genome_dp20`,
        `q`.`pct_genome_dp30` AS `pct_genome_dp30`,
        `q`.`vmr_depth` AS `vmr_depth`,
        `q`.`depth_q10` AS `depth_q10`,
        `q`.`depth_q20` AS `depth_q20`,
        `q`.`depth_q30` AS `depth_q30`,
        `q`.`raw_base_gb` AS `raw_base_gb`,
        `q`.`pct_overlap_reads` AS `pct_overlap_reads`,
        `q`.`pct_overlap_bases` AS `pct_overlap_bases`,
        `q`.`isize_iqr` AS `isize_iqr`,
        `q`.`isize_stdev` AS `isize_stdev`,
        `q`.`gc_depth_0_1` AS `gc_depth_0_1`,
        `q`.`gc_depth_1_5` AS `gc_depth_1_5`,
        `q`.`gc_depth_5_25` AS `gc_depth_5_25`,
        `q`.`gc_depth_25_75` AS `gc_depth_25_75`,
        `q`.`gc_depth_75_95` AS `gc_depth_75_95`,
        `q`.`gc_depth_95_99` AS `gc_depth_95_99`,
        `q`.`gc_depth_99_100` AS `gc_depth_99_100`,
        `q`.`library_size_m` AS `library_size_m`,
        
        `q`.`qc_pass` AS `qc_pass`,             // New rules from Tom, Nov 13, 2017
        `q`.`qc_flagged` AS `qc_flagged`,
        `q`.`qc_fail` AS `qc_fail`,
        
        FROM_UNIXTIME(`b`.`datearrived`) AS `recieved`,
        `b`.`bamsize` AS `size`,
        `s_qplot`.`name` AS `status_qplot`,
        `s_b37`.`name` AS `status_remap_hg37`,
        `s_b38`.`name` AS `status_remap_hg38`,
        `b`.`datemapping_b37` AS `mapped_b37`,
        `b`.`datemapping_b38` AS `mapped_b38`
    FROM
        ((((((`bamfiles` `b`
        JOIN `runs` `r` ON ((`b`.`runid` = `r`.`runid`)))
        JOIN `centers` `c` ON ((`r`.`centerid` = `c`.`centerid`)))
        JOIN `states` `s_qplot` ON ((`b`.`state_qplot` = `s_qplot`.`id`)))
        LEFT JOIN `qc_results` `q` ON ((`q`.`bam_id` = `b`.`bamid`)))
        JOIN `states` `s_b37` ON ((`b`.`state_b37` = `s_b37`.`id`)))
        JOIN `states` `s_b38` ON ((`b`.`state_b38` = `s_b38`.`id`)));

        0 AS `qc`,                              // Start of removed code
        IF(((`q`.`pct_freemix` < 3)
                AND (`q`.`pct_genome_dp10` > 95)
                AND (`q`.`mean_depth` > 30)),
            1,
            0) AS `qc_pass`,
        IF((`q`.`mean_depth` < 30), 1, 0) AS `qc_flagged`,
        IF(((`q`.`pct_freemix` > 3)
                OR (`q`.`pct_genome_dp10` < 95)),
            1,
            0) AS `qc_fail`,
                                                // End of removed code
*/

