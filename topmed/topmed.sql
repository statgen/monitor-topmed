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
mysql -h f-db -p    (Use FuS password)
flush hosts;
#   Verify it works by connecting on host that is having the problem

#   Clone the old database to the new
mysqldump -u sqlnhlbi --password=g9X+6iaO -h f-db nhlbi > /tmp/nhlbi.sql
mysql -u sqlnhlbi --password=g9X+6iaO -h localhost nhlbi < /tmp/nhlbi.sql
*/

/* --------------------------------------------------------------
#   Do not delete columns without telling Chris so he can fix mapping
#--------------------------------------------------------------- */
/* This table is used to control when/if an operation is permitted
   If an entry is in this table, it means it is disabled */
DROP TABLE IF EXISTS permissions;
CREATE TABLE permissions (
  id           INT         NOT NULL AUTO_INCREMENT,
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
  centerid     INT         NOT NULL AUTO_INCREMENT,
  centername   VARCHAR(16) NOT NULL,
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
  runid        INT         NOT NULL AUTO_INCREMENT,
  centerid     INT,
  dirname      VARCHAR(64) NOT NULL,
  status       VARCHAR(256),
  bamcount     INT,
  datayear     INT DEFAULT 2,          /* Year of project: 1, 2 ... */
  offsite      CHAR(1) DEFAULT('N'),   /* Original files kept offsite (N,Y,D) */
  xmlfound     INT DEFAULT 0,          /* Remove this */
  arrived      CHAR(1) DEFAULT('N'),   /* Y or N that all files arrived for this run */

  dateinit     VARCHAR(12),
  datecomplete VARCHAR(12),

  comments     TEXT,
  PRIMARY KEY  (runid)
);
CREATE INDEX index_centerid_dirname ON runs(centerid,dirname);
CREATE INDEX index_dirname ON runs(dirname);
 
/* Lists all BAMs we know about (e.g. NWDID)  Must be tied to a run */
DROP TABLE IF EXISTS bamfiles;
CREATE TABLE bamfiles (
  dateinit     VARCHAR(12),
  datayear     INT DEFAULT 2,           /* Year of project: 1, 2 ... */
  build        VARCHAR(4) DEFAULT '37', /* Build original input file user, 37, 38 etc */
  bamid        INT         NOT NULL AUTO_INCREMENT,
  runid        INT         NOT NULL,
  bamname_orig VARCHAR(96) NOT NULL,
  bamname      VARCHAR(96) NOT NULL,
  bamsize      VARCHAR(16) DEFAULT 0,
  nominal_length  INT DEFAULT 0,
  nominal_sdev INT DEFAULT 0,
  base_coord   INT DEFAULT 0,
  library_name VARCHAR(96),
  cramname     VARCHAR(96) NOT NULL,
  cramchecksum VARCHAR(96) NOT NULL,
  b37cramchecksum VARCHAR(96) NOT NULL,
  b38cramchecksum VARCHAR(96) NOT NULL,
  bamflagstat  BIGINT UNSIGNED DEFAULT NULL,
  cramflagstat BIGINT UNSIGNED DEFAULT NULL,
  b37flagstat  BIGINT UNSIGNED DEFAULT NULL,
  b38flagstat  BIGINT UNSIGNED DEFAULT NULL,
  studyname    VARCHAR(96) NOT NULL,
  piname       VARCHAR(96),
  phs          VARCHAR(12),
  phs_consent_short_name VARCHAR(24),
  phs_sra_sample_id VARCHAR(24),
  phs_sra_data_details VARCHAR(255),
  ncbierr      VARCHAR(255),
  checksum     VARCHAR(96) NOT NULL,
  expt_sampleid VARCHAR(24) DEFAULT 'UNKNOWN',  /* NWDID */
  nwdid_known  CHAR(1) DEFAULT 'N',

/* Fields to track state for each step */
/*
my $NOTSET    = 0;            # Not set
my $REQUESTED = 1;            # Task requested
my $SUBMITTED = 2;            # Task submitted to be run
my $STARTED   = 3;            # Task started
my $DELIVERED = 19;           # Data delivered, but not confirmed
my $COMPLETED = 20;           # Task completed successfully
my $CANCELLED = 89;           # Task cancelled
my $FAILEDCHECKSUM = 98;      # Task failed, because checksum at NCBI bad
my $FAILED    = 99;           # Task failed
*/ 
  state_arrive   INT DEFAULT 0,
  state_md5ver   INT DEFAULT 0,
  state_backup   INT DEFAULT 0,
  state_cram     INT DEFAULT 0,
  state_bai      INT DEFAULT 0,
  state_qplot    INT DEFAULT 0,
  state_b37      INT DEFAULT 0,
  state_b38      INT DEFAULT 0,

  state_ncbiexpt INT DEFAULT 0,     /* Experiment defined at NCBI (X) */
  state_ncbiorig INT DEFAULT 0,     /* Original input file as bam or cram (S) */
  state_ncbib37  INT DEFAULT 0,     /* Remapped cram build 37 as cram (P) */
  state_ncbib38  INT DEFAULT 0,     /* Remapped cram build 38 as cram (T) */
  time_ncbiexpt  CHAR(19),          /* yyy-mm-ssThh:hh:ss when loaded */
  time_ncbiorig  CHAR(19),
  time_ncbib37   CHAR(19),
  time_ncbib38   CHAR(19),

  state_gce38push  INT DEFAULT 0,   /* Handling b38 data in Google Cloud */
  state_gce38pull  INT DEFAULT 0,
  state_gce38post  INT DEFAULT 0,   /* Post process state_gce38pull data */

  PRIMARY KEY  (bamid)

/*   Handy queries

select bamid,bamname,runid,(select dirname from runs where runid=bamfiles.runid) from bamfiles where state_ncbib37=19;


  How can we have dozens of files with this checksum ?
select bamid,bamname from bamfiles where cramb37checksum='d41d8cd98f00b204e9800998ecf8427e'
*/
 
/* ####################################################
   Remove these some day 
   #################################################### */
  datearrived  VARCHAR(12),
  datemd5ver   VARCHAR(12),
  datebackup   VARCHAR(12),
  datecram     VARCHAR(12),
  datebai      VARCHAR(12),
  dateqplot    VARCHAR(12),
  datecp2ncbi  VARCHAR(12),
  datemapping  VARCHAR(12),
  jobidarrived VARCHAR(12),
  jobidmd5ver  VARCHAR(12),
  jobidbackup  VARCHAR(12),
  jobidcram    VARCHAR(12),
  jobidbai     VARCHAR(12),
  jobidqplot   VARCHAR(12),
  jobidnwdid   VARCHAR(12),
  jobidncbiorig VARCHAR(12),
  jobidncbib37  VARCHAR(12),
  jobidncbib38  VARCHAR(12),
  jobidb37     VARCHAR(12),
  jobidb38     VARCHAR(12),
  studyid      INT         NOT NULL,
  cramorigsent CHAR(1) DEFAULT 'N',
#   Added for remapping with build37
  cramb37sent CHAR(1) DEFAULT 'N',
  cramb37checksum VARCHAR(96) NOT NULL,
  cramb38checksum VARCHAR(96) NOT NULL,
  bam_delivered VARCHAR(12),
  jobidcp2ncbi VARCHAR(12),
  jobidmapping VARCHAR(12),

  refname      VARCHAR(96) DEFAULT 'UNKNOWN',
  expt_refname VARCHAR(96) DEFAULT 'UNKNOWN',

);
CREATE INDEX index_runid   ON bamfiles(runid);
CREATE INDEX index_nwdid   ON bamfiles(expt_sampleid);
CREATE INDEX index_refname ON bamfiles(refname);

/* ALTER TABLE bamfiles ADD COLUMN datebai VARCHAR(12) AFTER datebackup; */
/* ALTER TABLE runs ADD  COLUMN offsitebackup CHAR(1) DEFAULT 'N' AFTER datayear */

/*  This table is keeps track of freezes - sets of samples 
*/
DROP TABLE IF EXISTS freezes;
CREATE TABLE freezes (
  id      INT PRIMARY KEY AUTO_INCREMENT,
  bamid   INT         NOT NULL,
  freeze  VARCHAR(45) NOT NULL,
  created_at   DATETIME  NOT NULL,  /* Cannot use DEFAULT CURRENT_TIMESTAMP, use plugin */
  modified_at  TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
  CONSTRAINT fk_freezes_1
    FOREIGN KEY (bamid)
    REFERENCES bamfiles (bamid)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION
);
CREATE INDEX index_bamid   ON freezes(bamid);
CREATE INDEX index_freeze  ON freezes(freeze);


/*  This table is regularly loaded from a tab delimited summary file from NCBI
    The NCBI data is loaded into this table to make it easier for us to garner
    the state of data sent there.

    mysql -h f-db.sph.umich.edu -u sqlnhlbi --password=g9X+6iaO --local-infile=1 nhlbi    
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
  count_b38         INT DEFAULT 0,
  avetime_b38       INT DEFAULT 0,
  ncbicount_b38     INT DEFAULT 0,

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
   View of QC _Metric data    nhlbi_qc_metrics
+-------------------+-------------+------+-----+---------------------+-------+
| Field             | Type        | Null | Key | Default             | Extra |
+-------------------+-------------+------+-----+---------------------+-------+
| id                | int(11)     | YES  |     | 0                   |       |
| bam_id            | int(11)     | YES  |     | NULL                |       |
| pct_freemix       | float       | YES  |     | NULL                |       |
| n_reads_m         | float       | YES  |     | NULL                |       |
| pct_mapped        | float       | YES  |     | NULL                |       |
| pct_mq0           | float       | YES  |     | NULL                |       |
| pct_paired        | float       | YES  |     | NULL                |       |
| pct_prop_paired   | float       | YES  |     | NULL                |       |
| mapped_gb         | float       | YES  |     | NULL                |       |
| q20_gb            | float       | YES  |     | NULL                |       |
| pct_q20_base      | float       | YES  |     | NULL                |       |
| mean_depth        | float       | YES  |     | NULL                |       |
| pct_genome_cov    | float       | YES  |     | NULL                |       |
| isize_mode        | float       | YES  |     | NULL                |       |
| isize_median      | float       | YES  |     | NULL                |       |
| pct_dups          | float       | YES  |     | NULL                |       |
| pct_genome_dp5    | float       | YES  |     | NULL                |       |
| pct_genome_dp10   | float       | YES  |     | NULL                |       |
| pct_genome_dp20   | float       | YES  |     | NULL                |       |
| pct_genome_dp30   | float       | YES  |     | NULL                |       |
| vmr_depth         | float       | YES  |     | NULL                |       |
| depth_q10         | float       | YES  |     | NULL                |       |
| depth_q20         | float       | YES  |     | NULL                |       |
| depth_q30         | float       | YES  |     | NULL                |       |
| raw_base_gb       | float       | YES  |     | NULL                |       |
| pct_overlap_reads | float       | YES  |     | NULL                |       |
| pct_overlap_bases | float       | YES  |     | NULL                |       |
| isize_iqr         | float       | YES  |     | NULL                |       |
| isize_stdev       | float       | YES  |     | NULL                |       |
| gc_depth_0_1      | float       | YES  |     | NULL                |       |
| gc_depth_1_5      | float       | YES  |     | NULL                |       |
| gc_depth_5_25     | float       | YES  |     | NULL                |       |
| gc_depth_25_75    | float       | YES  |     | NULL                |       |
| gc_depth_75_95    | float       | YES  |     | NULL                |       |
| gc_depth_95_99    | float       | YES  |     | NULL                |       |
| gc_depth_99_100   | float       | YES  |     | NULL                |       |
| library_size_m    | float       | YES  |     | NULL                |       |
| created_at        | datetime    | YES  |     | NULL                |       |
| modified_at       | timestamp   | YES  |     | 0000-00-00 00:00:00 |       |
| sample_id         | varchar(96) | YES  |     | UNKNOWN             |       |
| study             | varchar(96) | NO   |     | NULL                |       |
| center            | varchar(16) | NO   |     | NULL                |       |
| recieved          | datetime    | YES  |     | NULL                |       |
| size              | varchar(16) | YES  |     | 0                   |       |
| status_qplot      | varchar(45) | NO   |     | NULL                |       |
| status_remap_hg37 | varchar(45) | NO   |     | NULL                |       |
| status_remap_hg38 | varchar(45) | NO   |     | NULL                |       |
| status_backup     | varchar(45) | NO   |     | NULL                |       |
| status_ncbi_hg37  | varchar(45) | NO   |     | NULL                |       |
| status_ncbi_hg38  | varchar(45) | NO   |     | NULL                |       |
+-------------------+-------------+------+-----+---------------------+-------+
   #################################################### */

/* ####################################################
   QC _Metric data    qc_results
+-------------------+-----------+------+-----+-------------------+----------------+
| Field             | Type      | Null | Key | Default           | Extra          |
+-------------------+-----------+------+-----+-------------------+----------------+
| id                | int(11)   | NO   | PRI | NULL              | auto_increment |
| bam_id            | int(11)   | NO   | MUL | NULL              |                |
| pct_freemix       | float     | YES  |     | NULL              |                |
| n_reads_m         | float     | YES  |     | NULL              |                |
| pct_mapped        | float     | YES  |     | NULL              |                |
| pct_mq0           | float     | YES  |     | NULL              |                |
| pct_paired        | float     | YES  |     | NULL              |                |
| pct_prop_paired   | float     | YES  |     | NULL              |                |
| mapped_gb         | float     | YES  |     | NULL              |                |
| q20_gb            | float     | YES  |     | NULL              |                |
| pct_q20_base      | float     | YES  |     | NULL              |                |
| mean_depth        | float     | YES  |     | NULL              |                |
| pct_genome_cov    | float     | YES  |     | NULL              |                |
| isize_mode        | float     | YES  |     | NULL              |                |
| isize_median      | float     | YES  |     | NULL              |                |
| pct_dups          | float     | YES  |     | NULL              |                |
| pct_genome_dp5    | float     | YES  |     | NULL              |                |
| pct_genome_dp10   | float     | YES  |     | NULL              |                |
| pct_genome_dp20   | float     | YES  |     | NULL              |                |
| pct_genome_dp30   | float     | YES  |     | NULL              |                |
| vmr_depth         | float     | YES  |     | NULL              |                |
| depth_q10         | float     | YES  |     | NULL              |                |
| depth_q20         | float     | YES  |     | NULL              |                |
| depth_q30         | float     | YES  |     | NULL              |                |
| raw_base_gb       | float     | YES  |     | NULL              |                |
| pct_overlap_reads | float     | YES  |     | NULL              |                |
| pct_overlap_bases | float     | YES  |     | NULL              |                |
| isize_iqr         | float     | YES  |     | NULL              |                |
| isize_stdev       | float     | YES  |     | NULL              |                |
| gc_depth_0_1      | float     | YES  |     | NULL              |                |
| gc_depth_1_5      | float     | YES  |     | NULL              |                |
| gc_depth_5_25     | float     | YES  |     | NULL              |                |
| gc_depth_25_75    | float     | YES  |     | NULL              |                |
| gc_depth_75_95    | float     | YES  |     | NULL              |                |
| gc_depth_95_99    | float     | YES  |     | NULL              |                |
| gc_depth_99_100   | float     | YES  |     | NULL              |                |
| library_size_m    | float     | YES  |     | NULL              |                |
| created_at        | datetime  | NO   |     | NULL              |                |
| modified_at       | timestamp | NO   |     | CURRENT_TIMESTAMP |                |
+-------------------+-----------+------+-----+-------------------+----------------+


/* ####################################################
   Remove these some day 
   #################################################### */
DROP TABLE IF EXISTS studies;
CREATE TABLE studies (
  studyid      INT         NOT NULL AUTO_INCREMENT,
  studyname    VARCHAR(64) NOT NULL,
  PRIMARY KEY  (studyid)
);
CREATE INDEX index_studyid_studyname ON studies(studyid,studyname);
CREATE INDEX index_studyname ON studies(studyname);

DROP TABLE IF EXISTS requestfiles;
CREATE TABLE requestfiles (
  reqid        INT         NOT NULL AUTO_INCREMENT,
  centerid     INT         NOT NULL,
  hostname     VARCHAR(255) NOT NULL,
  fetchpath    VARCHAR(255) NOT NULL,
  daterequestor VARCHAR(12),
  iprequestor  VARCHAR(16),
  ipuser       VARCHAR(16),
  status       VARCHAR(127),
  PRIMARY KEY  (reqid)
);
CREATE INDEX index_reqid ON requestfiles(centerid);


/* Table of data received from Adam Stine - I think this shows all files they have */
DROP TABLE IF EXISTS adam;
CREATE TABLE adam (
  id           INT         NOT NULL AUTO_INCREMENT,
  name_1       VARCHAR(64),
  md5sum_1     VARCHAR(32),
  upload_id_1  INT,
  track_date_1 VARCHAR(8),
  status1      VARCHAR(16),
  name_2       VARCHAR(64),
  md5sum_2     VARCHAR(32),
  upload_id_2  INT,
  track_date_2 VARCHAR(8),
  status_2      VARCHAR(16),
  PRIMARY KEY  (id)
);


