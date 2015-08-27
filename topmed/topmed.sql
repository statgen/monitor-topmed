#---------------------------------------------------------------
#       NHLBI TopMed tracking entries
#---------------------------------------------------------------
DROP TABLE IF EXISTS centers;
CREATE TABLE centers (
  centerid     INT         NOT NULL AUTO_INCREMENT,
  centername   VARCHAR(16) NOT NULL,
  PRIMARY KEY  (centerid)
);
INSERT INTO centers (centername) VALUES('broad');
INSERT INTO centers (centername) VALUES('illumina');
INSERT INTO centers (centername) VALUES('nygc');
INSERT INTO centers (centername) VALUES('uw');

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
  xmlfound     INT DEFAULT 0,

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
  studyid      INT         NOT NULL,
  bamname      VARCHAR(96) NOT NULL,
  studyname    VARCHAR(96) NOT NULL,
  piname       VARCHAR(96),
  nwdid        VARCHAR(32),
  phs          VARCHAR(12),
  checksum     VARCHAR(96) NOT NULL,
  refname      VARCHAR(96) NOT NULL,
  expt_refname VARCHAR(96) NOT NULL,
  expt_sampleid VARCHAR(24) NOT NULL,
  datearrived  VARCHAR(12),
  datemd5ver   VARCHAR(12),
  datebackup   VARCHAR(12),
  datecram     VARCHAR(12),
  datebai      VARCHAR(12),
  dateqplot    VARCHAR(12),
  datemapping  VARCHAR(12),
  datecp2ncbi  VARCHAR(12),
  jobidarrived VARCHAR(12),
  jobidmd5ver  VARCHAR(12),
  jobidbackup  VARCHAR(12),
  jobidcram    VARCHAR(12),
  jobidbai     VARCHAR(12),
  jobidqplot   VARCHAR(12),
  jobidmapping VARCHAR(12),
  jobidcp2ncbi VARCHAR(12),
  bamsize      VARCHAR(16) DEFAULT 0,

  dateinit     VARCHAR(12),
  PRIMARY KEY  (bamid)
);
CREATE INDEX index_runid   ON bamfiles(runid);
CREATE INDEX index_refname ON bamfiles(refname);

# ALTER TABLE bamfiles ADD COLUMN datebai VARCHAR(12) AFTER datebackup;

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


# Remove everything in database for testing
#delete from centers;
#delete from runs;
#delete from studies;
#delete from bamfiles;
#delete from requestfiles;
#delete from reportdata;
#alter table bamfiles add column datereport  VARCHAR(12) after dateqplot

#   Irina chose to generate report data with her own code and process   2015 June
#DROP TABLE IF EXISTS reportdata;
#CREATE TABLE reportdata (
#  reportid     INT         NOT NULL AUTO_INCREMENT,
#  bamid        INT         NOT NULL,
#  vcfname      VARCHAR(96),
#  bamfile_remapped    VARCHAR(96),
#  vcffile_harmonized  VARCHAR(96),
#  qplotfile    VARCHAR(96),
#  totalreads   INT,
#  mappingrate  FLOAT,
#  maprate_mqpass FLOAT,
#  targetmapping  FLOAT,
#  zeromapqual  FLOAT,
#  mapqual      FLOAT,
#  pairedreads  FLOAT,
#  properpair   FLOAT,
#  mappedbases  FLOAT,
#  q20bases     FLOAT,
#  q20basespct  FLOAT,
#  meandepth    FLOAT,
#  genomecover  FLOAT,
#  eps_mse      FLOAT,
#  eps_cycle_mean FLOAT,
#  gcbiasmse    FLOAT,
#  isize_mode   FLOAT,
#  isize_medium FLOAT,
#  duptate      FLOAT,
#  qcfailrate   FLOAT,
#  base_comp_a  FLOAT,
#  base_comp_c  FLOAT,
#  base_comp_g  FLOAT,
#  base_comp_t  FLOAT,
#  base_comp_o  FLOAT,
#  PRIMARY KEY  (reportid)
#);
#CREATE INDEX index_reportid ON reportdata(reportid);
#CREATE INDEX index_bamid ON reportdata(bamid,reportid);

