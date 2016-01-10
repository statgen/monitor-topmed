/* --------------------------------------------------------------
#       NHLBI TopMed tracking entries
#
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

/* Lists all the centers we get BAMs from */
DROP TABLE IF EXISTS centers;
CREATE TABLE centers (
  centerid     INT         NOT NULL AUTO_INCREMENT,
  centername   VARCHAR(16) NOT NULL,
  centerdesc   VARCHAR(96) NOT NULL,
  designdesc   TEXT NOT NULL,
  PRIMARY KEY  (centerid)
);
INSERT INTO centers (centername,centerdesc,designdesc) VALUES('broad', 'Broad Institute','Illumina sequencing of Homo sapiens via random selection');
INSERT INTO centers (centername,centerdesc,designdesc) VALUES('illumina', 'Illumina Fast Track Services','PCR-Free Paired-end libraries are manually generated from 500ng-1ug of gDNA using the Illumina TruSeq DNA Sample Preparation Kit (Catalog #: FC-121-2001), based on the protocol in the TruSeq DNA PCR-Free Sample Preparation Guide.  Pre-fragmentation gDNA cleanup is performed using paramagnetic sample purification beads (Agencourt (TM) AMPure (TM) XP reagents, Beckman Coulter).  Samples are fragmented and libraries are size selected following fragmentation and end-repair using paramagnetic sample purification beads, targeting short inserts.  Final libraries are quality controlled for size using a gel electrophoretic separation system and are quantified.');
INSERT INTO centers (centername,centerdesc,designdesc) VALUES('nygc', 'New York Genome Center','Whole genome sequencing using Illumina TruSeq PCR-free DNA library preparation with 500ng input DNA, sequenced to >30x mean coverage with 2x150bp reads on HiSeq X.');
INSERT INTO centers (centername,centerdesc,designdesc) VALUES('uw', 'University of Washington Genome Sciences','equivalent to Illumina TruSeq PCR-free DNA sample prep');

/* Lists all runs (directories of data). Must be tied to a center */
DROP TABLE IF EXISTS runs;
CREATE TABLE runs (
  runid        INT         NOT NULL AUTO_INCREMENT,
  centerid     INT,
  dirname      VARCHAR(64) NOT NULL,
  status       VARCHAR(256),
  bamcount     INT,
  xmlfound     INT DEFAULT 0,          /* Remove this */

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
  b37bamchecksum VARCHAR(96) NOT NULL,
  b38bamchecksum VARCHAR(96) NOT NULL,
  studyname    VARCHAR(96) NOT NULL,
  piname       VARCHAR(96),
  phs          VARCHAR(12),
  phs_consent_short_name VARCHAR(24),
  phs_sra_sample_id VARCHAR(24),
  phs_sra_data_details VARCHAR(255),
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
  state_ncbiexpt INT DEFAULT 0,
  state_ncbiorig INT DEFAULT 0,
  state_ncbib37  INT DEFAULT 0,
  state_ncbib38  INT DEFAULT 0,

  PRIMARY KEY  (bamid)


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

  bamcount          INT DEFAULT 0,      /* Count of arrived bams */
  errcount          INT DEFAULT 0,      /* Count of bams sent to NCBI in error */
  loadedbamcount    INT DEFAULT 0,      /* Count of loaded BAMs at NCBI */
  PRIMARY KEY  (yyyymmdd)
);




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
