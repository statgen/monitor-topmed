/* --------------------------------------------------------------
#       NHLBI TOPMed tracking entries for RNA data
#--------------------------------------------------------------- 
*/

/* Lists all directories of data (e.g. runs). Must be tied to a center */
DROP TABLE IF EXISTS tx_projects;
CREATE TABLE tx_projects (
   rnaprojectid  INT NOT NULL AUTO_INCREMENT,
   centerid     INT,
   dirname      VARCHAR(256) NOT NULL,   	/* Path to where files can be found */
   count        INT,                        /* Count of RNA samples in this project */
   datayear     INT DEFAULT 2019,           /* Year when data arrived, cannot use YEAR(NOW()) */
   arrived      CHAR(1) DEFAULT 'N',        /* Y or N that all files arrived */
   status       VARCHAR(256),
   dateinit     VARCHAR(12),
   datecomplete VARCHAR(12),
   comments     TEXT,

   PRIMARY KEY  (rnaprojectid)
);
/* Not sure these are needed */
/* CREATE INDEX index_rnaprojectid_rnaprojectpath ON tx_projects(rnaprojectid,projectpath); */
/* CREATE INDEX index_rnaprojectpath ON tx_projects(projectpath); */

/* All samples and their data. Must be tied to a tx_projects */
DROP TABLE IF EXISTS tx_samples;
CREATE TABLE tx_samples (
   txseqid       INT         NOT NULL AUTO_INCREMENT,
   rnaprojectid  INT         NOT NULL,      /* tx_projects id  */
   expt_sampleid VARCHAR(24),               /* NWDID or really TOR ID */
   fileprefix    VARCHAR(64),     		    /* Prefix of all filenames for this sample */
   count         INT  DEFAULT 0,            /* How many files for this sample */
   rnasubject    VARCHAR(64),               /* Investigator ID */
   notes         TEXT,                      /* Notes from tab/xls file */
   piname        VARCHAR(96),
   dateinit      VARCHAR(12),
   datearrived   VARCHAR(12),				/* Needed for tracking state of steps */
   emsg          VARCHAR(255),
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
   state_backup   INT DEFAULT 0,
   state_aws38copy INT DEFAULT 0,    /* Copy local data to AWS bucket */
   state_fix INT DEFAULT 0,          /* Track efforts to fix screwups */

   PRIMARY KEY  (txseqid)
);
CREATE INDEX index_rnatxseqid     ON tx_samples(rnaprojectid);
CREATE INDEX index_rnatx_sampleid ON tx_samples(expt_sampleid);

/* Each sample can have a number of associated files */
DROP TABLE IF EXISTS tx_files;
CREATE TABLE tx_files (
   fileid        INT         PRIMARY KEY AUTO_INCREMENT,
   txseqid       INT         NOT NULL,     /* Which tx_samples */
   filename      VARCHAR(256) NOT NULL,
   checksum      VARCHAR(32) DEFAULT ' ',
   intar         CHAR(1) DEFAULT 'N',      /* Is this file in tar file */
   dateinit      VARCHAR(12),

   CONSTRAINT unique_txseqid_filename UNIQUE (txseqid, filename)
);
CREATE INDEX index_rnatxfileid ON tx_files(fileid);
CREATE INDEX index_fileid_txseqid ON tx_files(fileid,txseqid);

/*   Handy queries
    ALTER TABLE tx_samples ADD COLUMN datebai VARCHAR(12) AFTER datebackup;
    ALTER TABLE tx_samples ADD  COLUMN offsitebackup CHAR(1) DEFAULT 'N' AFTER datayear;
    ALTER TABLE tx_samples MODIFY COLUMN datayear INT DEFAULT 5;
    ALTER TABLE tx_samples DROP COLUMN colname;
*/

