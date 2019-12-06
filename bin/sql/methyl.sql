/* --------------------------------------------------------------
#       NHLBI TOPMed tracking entries for Methylation data
#--------------------------------------------------------------- 
*/

/* Lists all directories of data (e.g. runs). Must be tied to a center */
DROP TABLE IF EXISTS methyl_projects;
CREATE TABLE methyl_projects (
   methylprojectid  INT NOT NULL AUTO_INCREMENT,
   centerid     INT,
   dirname      VARCHAR(256) NOT NULL,   	/* Path to where batch directies can be found */
   datayear     INT DEFAULT 2019,           /* Year when data arrived, cannot use YEAR(NOW()) */
   arrived      CHAR(1) DEFAULT 'N',        /* Y or N that all files arrived */
   status       VARCHAR(256),
   dateinit     VARCHAR(12),
   datecomplete VARCHAR(12),
   comments     TEXT,

   PRIMARY KEY  (methylprojectid)
);

/* All data files. Must be tied to a methyl_projects */
DROP TABLE IF EXISTS methyl_batch;
CREATE TABLE methyl_batch (
   methylbatchid    INT      NOT NULL AUTO_INCREMENT,
   methylprojectid  INT      NOT NULL,      /* Project from above */
   batchname     VARCHAR(24),               /* E.g. CRA-B0001 */
   count        INT,                        /* Count of Methylation samples in this batch */
   /* Sample is really four zip files */
   file0         VARCHAR(64),               /* Name of file */
   file0checksum VARCHAR(32) DEFAULT ' ',
   file0size     INT,
   file1         VARCHAR(64),
   file1checksum VARCHAR(32) DEFAULT ' ',
   file1size     INT,
   file2         VARCHAR(64),
   file2checksum VARCHAR(32) DEFAULT ' ',
   file2size     INT,
   file3         VARCHAR(64),
   file3checksum VARCHAR(32) DEFAULT ' ',
   file3size     INT,
   /* Columns to manage actions by the monitor */
   dateinit      VARCHAR(12),               /* When project data arrived */
   datearrived   VARCHAR(12),				/* Date of level0 file */
   emsg          VARCHAR(255),
   /* Fields to track state for each step. Uses same values as in bamfiles state_* */
   state_arrive   INT DEFAULT 0,
   state_verify   INT DEFAULT 0,
   state_backup   INT DEFAULT 0,
   state_aws38copy INT DEFAULT 0,    		/* Copy local data to AWS bucket */
   state_fix INT DEFAULT 0,          		/* Track efforts to fix screwups */

   PRIMARY KEY  (methylbatchid)
);
CREATE INDEX index_rnamethylbatchid ON methyl_batch(methylprojectid);
CREATE INDEX index_rnabatchname ON methyl_batch(batchname);

/* All samples must be tied to a methyl_batch */
DROP TABLE IF EXISTS methyl_samples;
CREATE TABLE methyl_samples (
   methylid   INT   NOT NULL AUTO_INCREMENT,
   methylbatchid    INT      NOT NULL,      /* methyl_batch id above */
   expt_sampleid VARCHAR(24),               /* NWDID or really TOE ID */
   array         CHAR(1),                   /* Array from Manifest */
   p05           FLOAT,                     /* p<0.05 from Manifest */
   p01           FLOAT,                     /* p<0.01 from Manifest */
   type          VARCHAR(8),                /* E.g. BIS or OXY from grn.idat file */
   version       VARCHAR(3),                /* E.g. v01 from grn.idat file */
   rowcol        CHAR(6),                   /* E.g. R03C01 from grn.idat file */
   methmatch     VARCHAR(24),               /* Future use */
   subject       VARCHAR(24),               /* Future use */
   /* Flags indicating if Manifest said various files existed
      Y - exists as indicated by single Manifest file
      C - exists as indicated by a combined Manifest file
      N - does not exist
   */
   grn_idat      CHAR(1) DEFAULT 'N',
   red_idat      CHAR(1) DEFAULT 'N',
   sdf           CHAR(1) DEFAULT 'N',
   meth_raw      CHAR(1) DEFAULT 'N',
   unmeth_raw    CHAR(1) DEFAULT 'N',
   beta_raw      CHAR(1) DEFAULT 'N',
   pvalue        CHAR(1) DEFAULT 'N',
   rgset_snp     CHAR(1) DEFAULT 'N',
   meth_noob     CHAR(1) DEFAULT 'N',
   unmeth_noob   CHAR(1) DEFAULT 'N',
   beta_noob     CHAR(1) DEFAULT 'N',
   noob_correct  CHAR(1) DEFAULT 'N',

   PRIMARY KEY  (methylid)
);
CREATE INDEX index_rnamethylid      ON methyl_samples(methylbatchid);
CREATE INDEX index_rnaexpt_sampleid ON methyl_samples(expt_sampleid);
CREATE UNIQUE INDEX index_methsampleid ON methyl_samples(expt_sampleid);
