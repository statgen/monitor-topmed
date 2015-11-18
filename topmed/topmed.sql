#---------------------------------------------------------------
#       NHLBI TopMed tracking entries
#
#   Do not delete columns without telling Chris so he can fix mapping
#---------------------------------------------------------------
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

#   The first of these,  EXPERIMENT -> alias,  identifies which record was pointed to by
#   RUN -> EXPERIMENT_REF -> refname. Then  EXPERIMENT -> SAMPLE_DESCRIPTOR -> refname
#   gives the individual identifier that we will search for in an external lookup table
#   to find out whether the associated .bam file should be transmitted immediately
#   to NCBI, or backed up for later transfer, or neither (if it already exists
#   in Amazon S3 or equivalent).
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

DROP TABLE IF EXISTS studies;
CREATE TABLE studies (
  studyid      INT         NOT NULL AUTO_INCREMENT,
  studyname    VARCHAR(64) NOT NULL,
  PRIMARY KEY  (studyid)
);
CREATE INDEX index_studyid_studyname ON studies(studyid,studyname);
CREATE INDEX index_studyname ON studies(studyname);

DROP TABLE IF EXISTS bamfiles;
CREATE TABLE bamfiles (
  bamid        INT         NOT NULL AUTO_INCREMENT,
  runid        INT         NOT NULL,
  studyid      INT         NOT NULL,          /* Remove this */
  bamname_orig VARCHAR(96) NOT NULL,
  bamname      VARCHAR(96) NOT NULL,
  bamsize      VARCHAR(16) DEFAULT 0,
  nominal_length  INT DEFAULT 0,
  nominal_sdev INT DEFAULT 0,
  base_coord   INT DEFAULT 0,
  library_name VARCHAR(96),
  cramname     VARCHAR(96) NOT NULL,
  cramchecksum VARCHAR(96) NOT NULL,
  studyname    VARCHAR(96) NOT NULL,
  piname       VARCHAR(96),
  phs          VARCHAR(12),
  phs_consent_short_name VARCHAR(24),
  phs_sra_sample_id VARCHAR(24),
  phs_sra_data_details VARCHAR(255),
  checksum     VARCHAR(96) NOT NULL,
  refname      VARCHAR(96) NOT NULL,
  expt_refname VARCHAR(96) NOT NULL,
  expt_sampleid VARCHAR(24) NOT NULL,
  nwdid_known  CHAR(1) DEFAULT 'N',
  datearrived  VARCHAR(12),
  datemd5ver   VARCHAR(12),
  datebackup   VARCHAR(12),
  datecram     VARCHAR(12),
  datebai      VARCHAR(12),
  dateqplot    VARCHAR(12),
  datecp2ncbi  VARCHAR(12),
  jobidarrived VARCHAR(12),
  jobidmd5ver  VARCHAR(12),
  jobidbackup  VARCHAR(12),
  jobidcram    VARCHAR(12),
  jobidbai     VARCHAR(12),
  jobidqplot   VARCHAR(12),
  jobidcp2ncbi VARCHAR(12),
  bam_delivered VARCHAR(12),
  cramorigsent CHAR(1) DEFAULT 'N',
#   Added for remapping with build37
  datemapping  VARCHAR(12),
  jobidmapping VARCHAR(12),
  cramb37sent CHAR(1) DEFAULT 'N',
  cramb37checksum VARCHAR(96) NOT NULL,
#   Added for remapping with build38
#  datemapping8 VARCHAR(12),
#  jobidmapping8 VARCHAR(12),
#  cramb39sent CHAR(1) DEFAULT 'N',
#  cramb39checksum VARCHAR(96) NOT NULL,

  dateinit     VARCHAR(12),
  PRIMARY KEY  (bamid)
);
CREATE INDEX index_runid   ON bamfiles(runid);
CREATE INDEX index_nwdid   ON bamfiles(expt_sampleid);
CREATE INDEX index_refname ON bamfiles(refname);

# ALTER TABLE bamfiles ADD COLUMN datebai VARCHAR(12) AFTER datebackup;


#   This table is used to control when/if an operation is permitted
#   If an entry is in this table, it means it is disabled
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

#   This table has never been used yet
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
