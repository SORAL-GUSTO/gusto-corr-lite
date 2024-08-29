#define _XOPEN_SOURCE
#include <inttypes.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include "callback.h"
#include "corrspec.h"
#include <errno.h>
#include <curl/curl.h>
#include "influx.h"

#include <Python.h>

#include <stdio.h>
#include <stdbool.h>
#include <fitsio.h>

#define PI 3.14159
#define BUFSIZE 256
#define DEBUG 0
#define DOQC 1

// Function to perform circular shift by -1  [UNUSED]
// functionally identical to the even/odd IQ[i+1] scheme
// when IQ[last] is set zero via Hann window
void circularShiftLeft(float arr[], int size) {
    if (size <= 1) return; // No shift needed for arrays of size 0 or 1

    float firstElement = arr[0]; // Store the first element

    // Shift all elements one position to the left
    for (int i = 0; i < size - 1; i++) {
        arr[i] = arr[i + 1];
    }

    // Place the first element at the end of the array
    arr[size - 1] = firstElement;
}


void get_git_commit_info(const char *filename, char *commit_info) {
    char command[BUFSIZ];
    FILE *fp;
    char hash[BUFSIZ];

    // Construct the command to get the commit hash
    snprintf(command, sizeof(command), "git log -1 --format=%%cd --date=format-local:'%%Y-%%m-%%d %%H:%%M:%%S %%Z' --pretty=format:\"Last commit %%h by %%an %%ad\" -- %s", filename);
    fp = popen(command, "r");
    if (fp == NULL) {
        perror("popen");
        exit(EXIT_FAILURE);
    }
    fgets(hash, sizeof(hash), fp);
    pclose(fp);
    hash[strcspn(hash, "\n")] = 0; // Remove trailing newline

    // Combine hash and date into the commit_info string
    snprintf(commit_info, BUFSIZ+sizeof(filename), "%s %s", filename, hash);
}

void get_proctime(char *proctime) {
   char command[BUFSIZ];
   FILE *fp;
   char date[BUFSIZ];
   
   snprintf(command, sizeof(command), "date +%%Y%%m%%d_%%H%%M%%S");
   fp = popen(command, "r");
   if (fp == NULL) {
      perror("popen");
      exit(EXIT_FAILURE);
   }
   fgets(date, sizeof(date), fp);
   pclose(fp);
   date[strcspn(date, "\n")] = 0; // Remove trailing newline

   snprintf(proctime, BUFSIZ, "%s", date);
}


void append_to_fits_table(const char *filename, struct s_header *fits_header, double *array) {
    fitsfile *fptr;  // FITS file pointer
		     //
    int status = 0;  // CFITSIO status value MUST be initialized to zero!
    int hdutype;
    int array_length;
    long nrows;
    char extname[] = "DATA_TABLE";
    long n_elements = sizeof(array) / sizeof(array[0]);

    CURL *curl;
    char *query = malloc(BUFSIZ);
    char *buf;

    char commit_info[256];
    char proctime[256];

    char line[4];
    float linefreq;
    int band, npix;

    // Try to open the FITS file in read/write mode. If it doesn't exist, create a new one.
    if (fits_open_file(&fptr, filename, READWRITE, &status)) {
        if (status == FILE_NOT_OPENED) {
            status = 0;  // Reset the status
            if (fits_create_file(&fptr, filename, &status)) {
                fits_report_error(stderr, status);  // Print any error message
                return;
            }

            // Create a primary array image (needed before any extensions can be created)
            if (fits_create_img(fptr, BYTE_IMG, 0, NULL, &status)) {
                fits_report_error(stderr, status);  // Print any error message
                return;
            }

            // Construct the primary FITS HEADER

	    // Various indices and keywords in the primary header that depend on Band #
            if (fits_header->unit == 6){ //ACS5 B1
	       strcpy(line, "NII");
	       band = 1;
	       npix = 512;
            }
            if (fits_header->unit == 4){ //ACS3 B2
	       strcpy(line, "CII");
	       band = 2;
	       npix = 1024;
            }

	    // Create some Primary header keyword value pairs and fill them from the current fits_header struct
            fits_write_key(fptr, TINT,      "CALID",   &fits_header->CALID,  "ID of correlator calibration", &status);
            fits_write_key(fptr, TSTRING, "TELESCOP",  "GUSTO",   "", &status);
            fits_write_key(fptr, TSTRING, "LINE",      &line,     "", &status);
            fits_write_key(fptr, TFLOAT,  "LINEFREQ",  &linefreq, "", &status);
            fits_write_key(fptr, TINT,    "BAND",      &band,     "GUSTO band #",     &status);
            fits_write_key(fptr, TINT,    "NPIX",      &npix,     "N spec pixels",    &status);
            fits_write_key(fptr, TSTRING, "DLEVEL",    "0.7",      "data level",      &status);
            get_proctime(proctime);
            fits_write_key(fptr, TSTRING, "Proctime",  proctime,  "processing time", &status);
            fits_write_key(fptr, TINT,    "SER_FLAG",  "0",       "SERIES FLAG",     &status);


            // Define the column parameters
            char *ttype[] ={"MIXER", "NINT", "UNIXTIME", "NBYTES", "CORRTIME", "INTTIME", "ROW_FLAG", "Ihigh", \
			    "Qhigh", "Ilow", "Qlow", "Ierr", "Qerr", "VIhi", "VQhi", "VIlo", "VQlo", "scanID", \
		            "subScan", "scan_type", "THOT", "RA", "DEC", "filename", "PSat", "Imon", "Gmon",   \
			    "spec", "CHANNEL_FLAG"};

            char *tunit[] ={" ", " ", "sec", " ", " ", "sec", " ", "counts", "counts", "counts",   \
			    "counts", " ", " ", "Volts", "Volts", "Volts", "Volts", " ", " ", " ", \
		            "Kelvin", "degrees", "degrees", "text", "Amps", "uA", "Amps", " ", " "};

	    // All header values are signed 32-bit except UNIXTIME which is uint64_t
            char *tform[29];
	    tform[0]  = "1J"; //int     mixer
	    tform[1]  = "1J"; //int	nint
	    tform[2]  = "1D"; //double	64 bit unixtime+fractional
	    tform[3]  = "1J"; //int	nbytes
	    tform[4]  = "1J"; //int     corrtime
	    tform[5]  = "1E"; //float   integration time
            tform[6]  = "1J"; //32 bit row flag
            tform[7]  = "1J"; //int	ihi
            tform[8]  = "1J"; //int	qhi
            tform[9]  = "1J"; //int	ilo
            tform[10] = "1J"; //int	qlo
            tform[11] = "1J"; //int	ierr
            tform[12] = "1J"; //int	qerr
            tform[13] = "1E"; //float	Vdac
            tform[14] = "1E"; //float	Vdac
            tform[15] = "1E"; //float	Vdac
            tform[16] = "1E"; //float	Vdac
            tform[17] = "1J"; //int	scanID
            tform[18] = "1J"; //int	subScan
            tform[19] = "6A"; //char    scan type
            tform[20] = "1E"; //float	THOT
            tform[21] = "1E"; //float	RA
            tform[22] = "1E"; //float	DEC
	    tform[23] = "48A";//char    filename
            tform[24] = "1E"; //float	PSat
            tform[25] = "1E"; //float	Imon
	    tform[26] = "1E"; //float   Gmon
	    // Various indices and keywords in the per row header that depend on Band #
            if (fits_header->unit==6){ //ACS5 B1
	       tform[27] = "512E";  //32 bit float
	       tform[28] = "512I";  //16 bit int
	    }
            if (fits_header->unit==4){ //ACS3 B2
	       tform[27] = "1024E";  //32 bit float
	       tform[28] = "1024I";  //16 bit int
	    }

	    int tfields = 29;

            // Create a binary table
            if (fits_create_tbl(fptr, BINARY_TBL, 0, tfields ,ttype, tform, tunit, extname, &status)) {
                fits_report_error(stderr, status);  // Print any error message
                return;
            }

        } else {
            fits_report_error(stderr, status);  // Print any error message
            return;
        }
    } ////////// END PRIMARY HEADER SECTION //////////

    // need to do this *again* because the one several lines back is only done for the fits file creation.
    // that bug took a day to find.
    if (fits_header->unit==6){ //ACS5 B1
        array_length=512;
    }
    if (fits_header->unit==4){ //ACS3 B2
        array_length=1024;
    }

    // Move to the named HDU (where the table is stored)
    if (fits_movnam_hdu(fptr, BINARY_TBL, extname, 0, &status)) {
        fits_report_error(stderr, status);  // Print any error message
        return;
    }

    // Get the current number of rows in the table
    if (fits_get_num_rows(fptr, &nrows, &status)) {
        fits_report_error(stderr, status);  // Print any error message
        return;
    }

    // insert a single empty row at the end of the output table
    if (fits_insert_rows(fptr, nrows, 1, &status)) {
        fits_report_error(stderr, status);  // Print any error message
        return;
    }

    // Write the header data
    fits_write_col(fptr, TINT32BIT,  1, nrows+1, 1, 1, &fits_header->mixer,    &status);
    fits_write_col(fptr, TINT32BIT,  2, nrows+1, 1, 1, &fits_header->nint,     &status);
    fits_write_col(fptr, TDOUBLE,    3, nrows+1, 1, 1, &fits_header->fulltime, &status);
    fits_write_col(fptr, TINT32BIT,  4, nrows+1, 1, 1, &fits_header->nbytes,   &status);
    fits_write_col(fptr, TINT32BIT,  5, nrows+1, 1, 1, &fits_header->corrtime, &status);
    fits_write_col(fptr, TFLOAT,     6, nrows+1, 1, 1, &fits_header->inttime,  &status);
    fits_write_col(fptr, TINT32BIT,  7, nrows+1, 1, 1, &fits_header->row_flag, &status);

    fits_write_col(fptr, TINT32BIT,  8, nrows+1, 1, 1, &fits_header->Ihi,      &status);
    fits_write_col(fptr, TINT32BIT,  9, nrows+1, 1, 1, &fits_header->Qhi,      &status);
    fits_write_col(fptr, TINT32BIT, 10, nrows+1, 1, 1, &fits_header->Ilo,      &status);
    fits_write_col(fptr, TINT32BIT, 11, nrows+1, 1, 1, &fits_header->Qlo,      &status);
    fits_write_col(fptr, TINT32BIT, 12, nrows+1, 1, 1, &fits_header->Ierr,     &status);
    fits_write_col(fptr, TINT32BIT, 13, nrows+1, 1, 1, &fits_header->Qerr,     &status);

    fits_write_col(fptr, TFLOAT,    14, nrows+1, 1, 1, &fits_header->VIhi,     &status);
    fits_write_col(fptr, TFLOAT,    15, nrows+1, 1, 1, &fits_header->VQhi,     &status);
    fits_write_col(fptr, TFLOAT,    16, nrows+1, 1, 1, &fits_header->VIlo,     &status);
    fits_write_col(fptr, TFLOAT,    17, nrows+1, 1, 1, &fits_header->VQlo,     &status);
									       
    fits_write_col(fptr, TINT32BIT, 18, nrows+1, 1, 1, &fits_header->scanID,   &status);
    fits_write_col(fptr, TINT32BIT, 19, nrows+1, 1, 1, &fits_header->subScan,  &status);
    fits_write_col(fptr, TSTRING,   20, nrows+1, 1, 1, &fits_header->type,     &status);

    fits_write_col(fptr, TSTRING,   24, nrows+1, 1, 1, &fits_header->filename, &status);

    // Write the spectra as a single 2*N column
    if (fits_write_col(fptr, TDOUBLE, 28, nrows+1, 1, 1 * array_length, array, &status)) {
        fits_report_error(stderr, status);  // Print any error message
        return;
    }
    // Write the channel flag as a single 2*N column
    if (fits_write_col(fptr, TDOUBLE, 29, nrows+1, 1, 1 * array_length, array, &status)) {
        fits_report_error(stderr, status);  // Print any error message
        return;
    }

    // Close the FITS file
    if (fits_close_file(fptr, &status)) {
        fits_report_error(stderr, status);  // Print any error message
        return;
    }

    if (DEBUG)
        printf("Array appended as a new row in the FITS table successfully.\n");
}

void printDateTimeFromEpoch(time_t ts)
{
   struct tm *tm = gmtime(&ts);

   char buffer[26];
   strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", tm);
   printf("UTC Date and Time: %s\n", buffer);
}


long get_seconds(const char *timestamp)
{
   struct tm t;
   char *pos;

   memset(&t, 0, sizeof(struct tm));

   strptime(timestamp, "%Y-%m-%dT%H:%M:%S", &t);

   time_t seconds = mktime(&t) - 25200;  //TODO: fix this hardcoded offset to UTC

   return (long) seconds;
}


int numPlaces (int n)
{
   if (n < 0)  return numPlaces ((n == INT_MIN) ? INT_MAX: -n);
   if (n < 10) return 1;
   return 1 + numPlaces (n / 10);
}


char nthdigit(int x, int n)
{
   while (n--) {
      x /= 10;
   }
   return (x % 10) + '0';
}


// Callback function to process the file
void const callback(char *filein){
   char *fullpath= malloc(48*sizeof(char));
   strcpy(fullpath, filein); // make a copy leaving filein intact for later tokenization

   char *datafile;	// datafile is filename with no path - used in fits header
   datafile = strrchr(fullpath, '/');
   if (datafile != NULL) {
	   datafile++;
   } else {
	   datafile = fullpath;
   }


   //char errfile[64] = "err.log";

   // InfluxDB easy_curl objects
   CURL *curl;
   CURLcode res;
   char *query = malloc(BUFSIZ);
   char *query_list[] = {"_VIhi", "_VQhi", "_VIlo", "_VQlo"};

   //timing
   struct timeval begin, end;
   gettimeofday(&begin, 0);

   //correlator file objectsI
   int N;
   FILE *fp;
   //FILE *fout;
   //FILE *errf;
   double P_I;
   double P_Q;
   double Ipwr;
   double Qpwr;
   struct corrType corr;

   // correlator variables from datafile
   uint64_t UNIXTIME;
   int NINT;
   int UNIT;
   int DEV;
   int NBYTES;
   float FRAC;
   int MIXER;

   float FS_FREQ; // Full scale frequency, B1==5000MHz || B2==5000MHz
   float dacV[4]; //0=VIhi, 1=VQhi, 2=VIlo, 3=VQlo
   float VIhi,  VQhi,  VIlo,  VQlo;

   // For integer correlator lags
   //int32_t *Rn;
   //int32_t *Rn2;
   
   // For normalized float correlator lags from Quant Corr
   float *Rn;
   float *Rn2;

   // file open notification
   printf("opened file: %s\n", filein);
   fp = fopen(filein, "r");

   // tokenize scanID from filename
   char *token;
   int position = 0;
   int band    = -1;
   char *prefix=malloc(7*sizeof(char));;
   int scanID  = -1;
   int subScan = -1;

   int CALID   = -1;	// scanID of previous correlator calibration
   bool error  = FALSE;	// status of error which may end processing early

   // Find file type from filename
   int i=0;
   char *ptr = NULL;
   const char *prefix_names[]={"SRC", "REF", "OTF", "HOT", "COLD", "FOC", "UNK"};

   // Use strtok to tokenize the filename using underscores as delimiters
   // First use "/" to redact any relative path
   //token = strtok(filein, "/");
   //token = strtok(NULL, "/");
   token = strtok(filein, "_");

   // Iterate through the tokens until reaching the 2nd position
   while (token != NULL ) {

      if (position == 0 ) {      //get band
	 if (strstr(token, "ACS3"))
            band = 2;
	 if (strstr(token, "ACS5"))
            band = 1;
         printf("Band is: %d\n", band);
      }

      if (position == 1 ) {      //get scan type
         while ( ptr==NULL ){
            ptr = strstr(token, prefix_names[i]);
            i++;
         }
         int len = strlen(prefix_names[i-1]);
         strncpy(prefix, ptr, len);
         prefix[len] = '\0';
         printf("The type is %s\n", prefix);
      }

      if (position == 2 ) {      //get scanID
         if (atoi(token)>0) scanID = atoi(token);
         printf("The scanID is: %d\n", scanID);
      }

      if (position == 3 ) {      //get subscanID
         subScan = atoi(token);
         printf("The data file index # is: %05d\n", subScan);
      }

      token = strtok(NULL, "_");
      position++;
   }

   // Build a regex with the range of the previous 16 scanID #s for Correlator DACs
   char scanIDregex[512];
   int pos = 0;
   pos += sprintf(&scanIDregex[pos], "^(");
   for (int k=0; k<15; k++){
      pos += sprintf(&scanIDregex[pos], "%d|", scanID-k);
   }
   pos += sprintf(&scanIDregex[pos], "%d)$", scanID-15);
   if (DEBUG)
      printf("%s\n", scanIDregex);

   // integer variables for fread from datafile
   int32_t value;
   uint32_t value1;	// for first 32 bits of UNIXTIMRE
   uint32_t value2;	//

   // figure out how many spectra in the file
   fseek(fp, 24, SEEK_SET);     // go to header position
   fread(&value, 4, 1, fp);
   int32_t bps = value;         // get bytes per spectra
   fseek(fp, 0L, SEEK_END);
   int sz = ftell(fp);          // get total bytes in file
   fseek(fp, 0L, SEEK_SET);     // go back to beginning of file
   printf("File has %.1f spectra\n", (float)sz/bps);

   int32_t header[22];
   const char *header_names[]={"UNIT", "DEV", "NINT", "UNIXTIME", "FRAC", "NBYTES", "CORRTIME", \
			       "EMPTY", "Ihigh", "Qhigh", "Ilow", "Qlow", "Ierr", "Qerr"};

//////////////////////////////  LOOP OVER ALL SPECTRA IN FILE  ///////////////////////////////////

   // Start at beginning of data file
   for (int j=0; j<(int)sz/bps; j++)
   {

   // Declare all Python objects for refpower
   // pArgs 
   PyObject *pArgs;
   PyObject *pValue;

   // Python objects for quantization correction
   PyObject *pArgsII, *pArgsQI, *pArgsIQ, *pArgsQQ; 
   PyObject *pListII, *pListQI, *pListIQ, *pListQQ; 
   PyObject *pValueII, *pValueQI, *pValueIQ, *pValueQQ; 

   // Arguments to relpower(XmonL, XmonH)
   if (DOQC){
      pArgs   = PyTuple_New(2); // for relpower(XmonL, XmonH)
      pArgsII = PyTuple_New(5); // for qc(ImonL, ImonH, ImonL, ImonH, corr.II)
      pArgsIQ = PyTuple_New(5); // for qc(ImonL, ImonH, QmonL, QmonH, corr.IQ)
      pArgsQI = PyTuple_New(5); // for qc(QmonL, QmonH, ImonL, ImonH, corr.IQ)
      pArgsQQ = PyTuple_New(5); // for qc(QmonL, QmonH, QmonL, QmonH, corr.QQ)
   }


   if (DEBUG)
      printf("The type is %s\n", prefix);
      // Loop over header location
      for (int i=0; i<22; i++){


         if (i==3){
            //UNIXTIME is 64 bits
            fread(&value1, 4, 1, fp); //Least significant 32 bits
            fread(&value2, 4, 1, fp); //Most significant 32 bis
            UNIXTIME = (((uint64_t)value2 << 32) | value1 ) / 1000.; //unixtime is to msec, store as 1sec int
            FRAC     = (((uint64_t)value2 << 32) | value1 ) % 1000 ; //fractional part is 1msec
         }
         else{
            fread(&value, 4, 1, fp);
         }
         header[i] = (value);
      }

      // fill variables from header array
      UNIT          = header[0];
      DEV           = header[1];
      NINT          = header[2];
      //UNIXTIME    = header[3]; 64 bits
      //CPU         = header[4]; broken?
      NBYTES        = header[5];
      corr.corrtime = header[6];
      //EMPTY         header[7];
      corr.Ihi      = header[8];
      corr.Qhi      = header[9];
      corr.Ilo      = header[10];
      corr.Qlo      = header[11];
      corr.Ierr     = header[12];
      corr.Qerr     = header[13];

      // Since we know UNIT and NLAGS, set correlator frequency and fits spectrum for later use
      // And also secret decoder ring (unit,dev) -> (Mixer #)
      int mixerTable[2][4] = {
	      {2, 3, 4, 6},	// band 1 mixers
	      {2, 3, 5, 8}	// band 2 mixers
      };
      MIXER = mixerTable[band-1][DEV-1];
      if (band==1){ //ACS5 B1
         FS_FREQ = 5000.;
      }
      if (band==2){ //ACS3 B2
         FS_FREQ = 5000.;
      }


      // Indication that this is a broken header file from correlator STOP signal
      if (corr.Ierr!=0 || corr.Qerr!=0 || \
           corr.Ihi==0 ||  corr.Qhi==0 || corr.Ilo==0 || corr.Qlo==0 || \
           (corr.corrtime*256.)/(FS_FREQ*1000000.)<0.1 || (corr.corrtime*256.)/(FS_FREQ*1000000.)>10.0 )
      {
         printf("######################## ERROR ###########################\n");
         printf("#                                                        #\n");
         printf("#                Error, data is no good!                 #\n");
         printf("#                        Exiting!                        #\n");
         printf("#                                                        #\n");
         printf("######################## ERROR ###########################\n");
         printf("CORRTIME was %.6f\n\n", (corr.corrtime*256.)/(FS_FREQ*1000000.));
	 error = TRUE;
         break;
      }

    
      /////////// Get data for the HIGH SPEED Housekeeping for per row fits header ///////////
      ///////////   This data *does* change every spectra   ///////////

      // Single SELECT for CORRELATOR DACS from current or nearest previous Correlator Cal
      curl = init_influx();
      sprintf(query, "&q=SELECT * FROM /^ACS%d_DEV%d_*/ WHERE \"scanID\"=~/^%s/ ORDER BY time DESC LIMIT 10", UNIT-1, DEV, scanIDregex);
      influxReturn = influxWorker(curl, query);
      dacV[0] = influxReturn->value[3]; //VIhi  Order is reverse alphabetical from InfluxDB SELECT
      dacV[1] = influxReturn->value[1]; //VQhi
      dacV[2] = influxReturn->value[2]; //VIlo
      dacV[3] = influxReturn->value[0]; //VQlo
      CALID = influxReturn->scanID;

      // Free Influx struct from ACS_DEV_VDAC
      freeinfluxStruct(influxReturn);

      // don't trust myself to rewrite the below, just copy from vector into floats
      VIhi = dacV[0];
      VQhi = dacV[1];
      VIlo = dacV[2];
      VQlo = dacV[3];

      // this section unfuck-ifys special cases when ICE was off by one
      if (VQlo==0.){
         VIhi=VIhi-(VIlo-VQhi);  //make up this lost data, it'l be close enough
         VQhi = dacV[0];
         VIlo = dacV[1];
         VQlo = dacV[2];
      }
      
      if (VIhi==0.){ //Still no values?  bail and don't make spectra
         printf("######################## ERROR ###########################\n");
         printf("#                                                        #\n");
         printf("#                  Error, no DAC values!                 #\n");
         printf("#                        Exiting!                        #\n");
         printf("#                                                        #\n");
         printf("######################## ERROR ###########################\n");
         break;
      }

      // DEBUG
      if (DEBUG)
         printf("VIhi %.3f\tVQhi %.3f\tVIlo %.3f\tVQlo %.3f\n", VIhi, VQhi, VIlo, VQlo);

      // Build the spectra filename and put it in the spectra directory
      //char fileout[512];

      //sprintf(fileout, "spectra/ACS%d_%s_%05d_DEV%d_INDX%04d_NINT%03d.txt", UNIT-1, prefix, scanID, DEV, subScan, j);
      //fout = fopen(fileout, "w");

      //read human readable "Number of Lags"
      if (NBYTES==8256)
         N = 512;
      else if (NBYTES==6208)
         N = 384;
      else if (NBYTES==4160)
         N = 256;
      else if (NBYTES==2112)
         N = 128;
      int specA = (int) N/128 - 1;

      //We don't know the lag # until we open the file, so malloc now
      corr.II   = malloc(N*sizeof(int32_t));   //Uncorrected ints
      corr.QI   = malloc(N*sizeof(int32_t));
      corr.IQ   = malloc(N*sizeof(int32_t));
      corr.QQ   = malloc(N*sizeof(int32_t));
      corr.IIqc = malloc(N*sizeof(float));     //Normalized Quantization Corrected floats
      corr.QIqc = malloc(N*sizeof(float));
      corr.IQqc = malloc(N*sizeof(float));
      corr.QQqc = malloc(N*sizeof(float));
      //Rn  = malloc(2*N*sizeof(int32_t));     //Rn,Rn2 ints for unquantization corrected
      //Rn2 = malloc(4*N*sizeof(int32_t));
      Rn  = malloc(2*N*sizeof(float));         //Rn,Rn2 floats for normalizaed, quantization corrected
      Rn2 = malloc(4*N*sizeof(float));

      // Read lags in from file in order after header
      for (int i=0; i<N; i++){
         fread(&value, 4, 1, fp);
         corr.II[i] = value;
         fread(&value, 4, 1, fp);
         corr.QI[i] = value;
         fread(&value, 4, 1, fp);
         corr.IQ[i] = value;
         fread(&value, 4, 1, fp);
         corr.QQ[i] = value;
      }



      // PASS THE UNCORRECTED LAGS INTO PYTHON FOR QUANTIZATION CORRECTION HERE !!!


      if (DOQC){
      // create a Python list and fill it with normalized uncorrected correlation coefficients
      // Check for errors
      if (pModule != NULL) {
         // Get the function from the module

         if (pFunc2 && PyCallable_Check(pFunc2)) {
            // For first four arguments: XmonL, XmonH, YmonL, YmonH
	    // II
	    PyTuple_SetItem(pArgsII, 0, PyFloat_FromDouble((double)corr.Ilo/(double)corr.corrtime));
	    PyTuple_SetItem(pArgsII, 1, PyFloat_FromDouble((double)corr.Ihi/(double)corr.corrtime));
	    PyTuple_SetItem(pArgsII, 2, PyFloat_FromDouble((double)corr.Ilo/(double)corr.corrtime));
	    PyTuple_SetItem(pArgsII, 3, PyFloat_FromDouble((double)corr.Ihi/(double)corr.corrtime));
            // IQ
	    PyTuple_SetItem(pArgsIQ, 0, PyFloat_FromDouble((double)corr.Ilo/(double)corr.corrtime));
	    PyTuple_SetItem(pArgsIQ, 1, PyFloat_FromDouble((double)corr.Ihi/(double)corr.corrtime));
	    PyTuple_SetItem(pArgsIQ, 2, PyFloat_FromDouble((double)corr.Qlo/(double)corr.corrtime));
	    PyTuple_SetItem(pArgsIQ, 3, PyFloat_FromDouble((double)corr.Qhi/(double)corr.corrtime));
            // QI
	    PyTuple_SetItem(pArgsQI, 0, PyFloat_FromDouble((double)corr.Qlo/(double)corr.corrtime));
	    PyTuple_SetItem(pArgsQI, 1, PyFloat_FromDouble((double)corr.Qhi/(double)corr.corrtime));
	    PyTuple_SetItem(pArgsQI, 2, PyFloat_FromDouble((double)corr.Ilo/(double)corr.corrtime));
	    PyTuple_SetItem(pArgsQI, 3, PyFloat_FromDouble((double)corr.Ihi/(double)corr.corrtime));
            // QQ
	    PyTuple_SetItem(pArgsQQ, 0, PyFloat_FromDouble((double)corr.Qlo/(double)corr.corrtime));
	    PyTuple_SetItem(pArgsQQ, 1, PyFloat_FromDouble((double)corr.Qhi/(double)corr.corrtime));
	    PyTuple_SetItem(pArgsQQ, 2, PyFloat_FromDouble((double)corr.Qlo/(double)corr.corrtime));
	    PyTuple_SetItem(pArgsQQ, 3, PyFloat_FromDouble((double)corr.Qhi/(double)corr.corrtime));

            // For fifth argument: convert C array to a Python List
	    // Creates new reference here.  pList** must be DECREF to free
            pListII = PyList_New(N);
            pListQI = PyList_New(N);
            pListIQ = PyList_New(N);
            pListQQ = PyList_New(N);
	   
            for (int i=0; i<N; ++i){
               PyList_SetItem(pListII, i, PyFloat_FromDouble((double)corr.II[i]/(double)corr.corrtime));
               PyList_SetItem(pListIQ, i, PyFloat_FromDouble((double)corr.IQ[i]/(double)corr.corrtime));
               PyList_SetItem(pListQI, i, PyFloat_FromDouble((double)corr.QI[i]/(double)corr.corrtime));
               PyList_SetItem(pListQQ, i, PyFloat_FromDouble((double)corr.QQ[i]/(double)corr.corrtime));
            }
	    // pList** are stolen by SetItem here.  Do not DECREF
	    PyTuple_SetItem(pArgsII, 4, pListII);
	    PyTuple_SetItem(pArgsIQ, 4, pListIQ);
	    PyTuple_SetItem(pArgsQI, 4, pListQI);
	    PyTuple_SetItem(pArgsQQ, 4, pListQQ);


	    // DEBUG print uncorrected
	    if (DEBUG) {
	       printf("XY   ");
               for (int i=0; i<10; ++i)
                  printf("%f ", (float) corr.II[i]/corr.II[0]);
	       printf("\n");
	    }
   
            // Call the function
            pValueII = PyObject_CallObject(pFunc2, pArgsII);
            pValueIQ = PyObject_CallObject(pFunc2, pArgsIQ);
            pValueQI = PyObject_CallObject(pFunc2, pArgsQI);
            pValueQQ = PyObject_CallObject(pFunc2, pArgsQQ);
   
            // Check for errors
            if (pValueII != NULL && pValueIQ != NULL && pValueQI != NULL && pValueQQ != NULL) {
               // Extract the values from the returned list
               for (int i=0; i<N; ++i){

		  // Fill corrqc struct with normalized float32 QC values
                  corr.IIqc[i] = PyFloat_AsDouble(PyList_GetItem(pValueII, i));
                  corr.IQqc[i] = PyFloat_AsDouble(PyList_GetItem(pValueIQ, i));
                  corr.QIqc[i] = PyFloat_AsDouble(PyList_GetItem(pValueQI, i));
                  corr.QQqc[i] = PyFloat_AsDouble(PyList_GetItem(pValueQQ, i));
               }
            } else{
	       PyErr_Print();
	    }


	    // DEBUG print corrected
	    if (DEBUG) {
	       printf("XYqc ");
               for (int i=0; i<10; ++i)
                  printf("%f ", corr.IIqc[i]);
	       printf("\n");
	    }


         } else {
	    if (PyErr_Occurred()) {
	       PyErr_Print();
	       fprintf(stderr, "Cannot find function\n");
	    }
	 }
      } else {
	 PyErr_Print();
	 fprintf(stderr, "Failed to load the Python module\n");
      }
   }  // END IF DOQC
      
      // If we're not going to do the QC DLL, then just copy lags over and normlize them at least
   else
   {
      for (int i=0; i<N; i++){
         corr.IIqc[i] = corr.II[i]/corr.II[0];
         corr.IQqc[i] = corr.IQ[i]/corr.IQ[0];
         corr.QIqc[i] = corr.QI[i]/corr.QI[0];
         corr.QQqc[i] = corr.QQ[i]/corr.QQ[0];
      }
   }  // END IF ELSE


      // GET THE QUANTIZATION CORRECTED LAGS BACK AND  CONTINUE ON WITH LAGS->SPECTRA


      // Combine IQ lags into R[]
         for (int i=0; i<(2*N)-1; i++){
            if (i%2 == 0) Rn[i] = 0.5* (corr.IIqc[i/2] + corr.QQqc[i/2]);           // Even indicies
            if (i%2 == 1) Rn[i] = 0.5* (corr.IQqc[(i-1)/2] + corr.QIqc[(i-1)/2+1]); // Odd  indicies
         }
         int le = 2*N-1; //last element
         Rn[le] = 0.5* (corr.IQqc[(le-1)/2] + corr.QIqc[(le-1)/2]);         // Last element

      // Mirror R[] symmetrically
         for (int i=0; i<4*N-1; i++){
           if ( i<2*N ) Rn2[2*N-1-i] = Rn[i];
           else        Rn2[i]       = Rn[i-(2*N-1)];
         }

/*
      // Try the Omnisys ICD spectrometer_api_CCv0005.docx routine
      //
      // Set QC[0] = 0
      corr.IQqc[0] = 0.;
      // circular shift minus one step
       circularShiftLeft(corr.IQqc, sizeof(corr.IQqc) / sizeof(corr.IQqc[0]));

        for (int i=0; i<2*N; i+=2){
           Rn[i]   = 0.5 * (corr.IIqc[i/2] + corr.QQqc[i/2]);
           Rn[i+1] = 0.5 * (corr.IQqc[i/2] + corr.QIqc[i/2]);
        }
      // For whatever reason, this doesn't work straight outta the box
*/

      // Apply Hann window
         for (int i=0; i<4*N; i++){
           Rn2[i] = Rn2[i]*0.5*(1-cos(2*PI*i/(4*N)));     //real with Hann window
         }

      // Fill fft struct
         for (int i=0; i<4*N; i++){
           spec[specA].in[i][0] = Rn2[i];  //real
           spec[specA].in[i][1] = 0.;      //imag
         }

      // Do FFT and print
         fftw_execute(spec[specA].p);

      // Compute Power Coefficients
         P_I = pow((VIhi-VIlo),2) * 1.8197 / pow((erfinv(1-2*(double)corr.Ihi/(double)corr.corrtime)) + \
                                                 (erfinv(1-2*(double)corr.Ilo/(double)corr.corrtime)),2);
         P_Q = pow((VQhi-VQlo),2) * 1.8197 / pow((erfinv(1-2*(double)corr.Qhi/(double)corr.corrtime)) + \
                                                 (erfinv(1-2*(double)corr.Qlo/(double)corr.corrtime)),2);

      // GET ALTERNATIVE POWER CALIBRATION FROM OMNISYS DLL
      // P_I,  P_Q calculated from the inverse error function using (Ihi, Ilow) is equal
      // to the relpower function from the Omnisys DLL to high precision. When  validated
      // over a range of power levels (HOT, OTF, REF) the two give equal results
      // (The factor 1.8197 likely comes from the erfinv for a 3-level ACS)
      
      if (DOQC){ // Only do this if Python Objects are declared
         PyTuple_SetItem(pArgs, 0, PyFloat_FromDouble((double)corr.Ilo/corr.corrtime));
         PyTuple_SetItem(pArgs, 1, PyFloat_FromDouble((double)corr.Ihi/corr.corrtime));
         pValue = PyObject_CallObject(pFunc1, pArgs);
         Ipwr = 2. * pow((VIhi-VIlo),2) * PyFloat_AsDouble(pValue);

         PyTuple_SetItem(pArgs, 0, PyFloat_FromDouble((double)corr.Qlo/corr.corrtime));
         PyTuple_SetItem(pArgs, 1, PyFloat_FromDouble((double)corr.Qhi/corr.corrtime));
         pValue = PyObject_CallObject(pFunc1, pArgs);
         Qpwr = 2. * pow((VQhi-VQlo),2) * PyFloat_AsDouble(pValue);


         // Output the relative power comparison
         if (DEBUG)
            printf("Ipwr is %f\tQpwr is %f\tP_I is %f\tP_Q is %f\n", Ipwr, Qpwr, P_I, P_Q);
      }

      // Comment out all text spectral file outputs.  FITS done below
      // Will be deprecated soon.
      // Header information in spectra file
/*
      fprintf(fout, "UNIXTIME\t%" PRIu64 "\n", UNIXTIME);
      fprintf(fout, "CORRTIME\t%.6f\n", (corr.corrtime*256.)/(FS_FREQ*1000000.));
      fprintf(fout, "UNIT\t%d\n", UNIT);
      fprintf(fout, "DEV\t%d\n",   DEV); 
      fprintf(fout, "NLAGS\t%d\n", N);
      fprintf(fout, "NINT\t%d\n", NINT);
      fprintf(fout, "LEVEL_Ihi\t%.2f\n",  (double)corr.Ihi/(double)corr.corrtime);
      fprintf(fout, "LEVEL_Qhi\t%.2f\n",  (double)corr.Qhi/(double)corr.corrtime);
      fprintf(fout, "LEVEL_Ilo\t%.2f\n",  (double)corr.Ilo/(double)corr.corrtime);
      fprintf(fout, "LEVEL_Qlo\t%.2f\n",  (double)corr.Qlo/(double)corr.corrtime);
      fprintf(fout, "IHI\t%u\n",  corr.Ihi);
      fprintf(fout, "QHI\t%u\n",  corr.Qhi);
      fprintf(fout, "ILO\t%u\n",  corr.Ilo);
      fprintf(fout, "QLO\t%u\n",  corr.Qlo);

    // Used calibration DAC voltages
      fprintf(fout, "VIHI\t%.3f\nVQHI\t%.3f\nVILO\t%.3f\nVQLO\t%.3f\n", VIhi, VQhi, VIlo, VQlo);
      fprintf(fout, "CAL\t%d\n", CALID);
      fprintf(fout, "DATA\t%d\n", scanID);

    // Data consistency checks and error flaging
      fprintf(fout, "IERR\t%u\n", corr.Ierr);
      fprintf(fout, "QERR\t%u\n", corr.Qerr);
      fprintf(fout, "ZEROLAGSUM_I\t%u\t%u\n", corr.Ihi+corr.Ilo, corr.II[0]);
      fprintf(fout, "ZEROLAGSUM_Q\t%u\t%u\n", corr.Qhi+corr.Qlo, corr.QQ[0]);
*/

      // Write to error log.
      // Was for flight ops.  Probably will deprecate soon.
      /*
      if (1)
      {
         errf = fopen(errfile, "a");
         fprintf(errf, "%s\t", name);
         fflush(errf);
         fprintf(errf, "%" PRIu64 "\t", UNIXTIME);
         fprintf(errf, "%u\t", CPU); //CPU => FRAC
         fprintf(errf, "%s\t", prefix);
         fprintf(errf, "%d\t%d\t", UNIT, DEV);
         fprintf(errf, "%.6f\t", (corr.corrtime*256.)/(FS_FREQ*1000000.));
         fprintf(errf, "%u\t", corr.Ierr);
         fprintf(errf, "%u\t", corr.Qerr);
         fprintf(errf, "%u\t", abs((corr.Ihi+corr.Ilo)-corr.II[0]));
         fprintf(errf, "%u\n", abs((corr.Qhi+corr.Qlo)-corr.QQ[0]));
         fclose(errf);
      }
      */

      //fprintf(fout, "ETAQ\t%.3f\n", 1/sqrt(P_I*P_Q));
   
      //Print in counts per second (assuming FS_FREQ sampling freq)
      /*
      for (int i=0; i<2*N; i++){
         fprintf(fout, "%d\t%lf\n", (FS_FREQ*i)/(2*N), sqrt(P_I*P_Q) * \
                                    sqrt(pow(spec[specA].out[i][0],2)+pow(spec[specA].out[i][1],2)));
      }
      fclose(fout); //close single spectra file
      */
		    

/*
 * Don't use this stuff, it was for V&V with the lags and spectra from LabView.
 * Never gets utilized now, will be deprecated in the future
      if (DEBUG>1){
         fout = fopen("lags.txt", "w");
         //Save QC lags
         for (int i=0; i<N; i++){
            fprintf(fout, "%d\t%f\t%f\t%f\t%f\n", i, corr.IIqc[i], corr.QQqc[i], corr.IQqc[i], corr.QIqc[i]);
         }
         fclose(fout);


         fout = fopen("lagsRn.txt", "w");
         //Save combined lags
         for (int i=0; i<2*N; i++){
            fprintf(fout, "%d\t%f\n", i,Rn[i]);
         }
         fclose(fout);


         fout = fopen("lagsRn2.txt", "w");
         //Save combined mirrored lags
         for (int i=0; i<4*N; i++){
            fprintf(fout, "%d\t%f\n", i,Rn2[i]);
         }
         fclose(fout);
      }
 */

      // Construct the per row FITS HEADER
      s_header *fits_header = (s_header *)malloc(sizeof(s_header));
      fits_header->type     = (char *)malloc(6);
      fits_header->filename = (char *)malloc(48);

      fits_header->unit     = UNIT;
      fits_header->dev      = DEV;
      fits_header->mixer    = MIXER;
      fits_header->nint     = NINT;
      fits_header->fulltime = (double) UNIXTIME + 0.001*FRAC;
      fits_header->nbytes   = NBYTES;
      fits_header->corrtime = (int) corr.corrtime;
      fits_header->inttime  =(corr.corrtime*256.)/(FS_FREQ*1000000.);
      fits_header->row_flag = 0;
      fits_header->channel_flag = 0;

      fits_header->Ihi      = corr.Ihi;
      fits_header->Qhi      = corr.Qhi;
      fits_header->Ilo      = corr.Ilo;
      fits_header->Qlo      = corr.Qlo;
      fits_header->Ierr     = corr.Ierr;
      fits_header->Qerr     = corr.Qerr;

      fits_header->VIhi     = VIhi;
      fits_header->VQhi     = VQhi;
      fits_header->VIlo     = VIlo;
      fits_header->VQlo     = VQlo;

      fits_header->scanID   = scanID;
      fits_header->subScan  = subScan;
      fits_header->CALID    = CALID;

      strcpy(fits_header->type, prefix);
      strcpy(fits_header->filename, datafile);
      
      // Construct the FITS DATA
      double *array = malloc(2*N*sizeof(double));
      for (int i=0; i<2*N; i++){
         array[i] = sqrt(P_I*P_Q) * sqrt(pow(spec[specA].out[i][0],2)+pow(spec[specA].out[i][1],2));
      }

      // Let's try out CFITSIO!
      char fitsfile[20];
      sprintf(fitsfile, "ACS%d_%s_%05d.fits", UNIT-1, prefix, scanID);
      if (DEBUG)
         printf("%s\n", fitsfile);
      append_to_fits_table(fitsfile, fits_header, array); 

      // Free items before next spectra within this file
      // All of these objects are malloced at the start of every spectra
      free(corr.II);      //free all mallocs
      free(corr.QI);
      free(corr.IQ);
      free(corr.QQ);
      free(corr.IIqc);    //free all mallocs
      free(corr.QIqc);
      free(corr.IQqc);
      free(corr.QQqc);
      free(Rn);
      free(Rn2);
      free(array);

      free(fits_header->type);
      free(fits_header->filename);
      free(fits_header);

      // Clean Up memory before moving to next spectra
      // All of these objects are malloced at the start of every spectra and freed every time
      // 4 X pArgs 
      //printf("free pArgsII\t from refcount\t %ld\n", pArgsII->ob_refcnt);
      Py_DECREF(pArgsII);
      //printf("free pArgsIQ\t from refcount\t %ld\n", pArgsIQ->ob_refcnt);
      Py_DECREF(pArgsIQ);
      //printf("free pArgsQI\t from refcount\t %ld\n", pArgsQI->ob_refcnt);
      Py_DECREF(pArgsQI);
      //printf("free pArgsQQ\t from refcount\t %ld\n", pArgsQQ->ob_refcnt);
      Py_DECREF(pArgsQQ);
    
      // and the refpower
      //printf("free pArgs\t from refcount\t %ld\n", pArgs->ob_refcnt);
      Py_DECREF(pArgs);

      // 4 X pValue
      //printf("free pValueII\t from refcount\t %ld\n", pValueII->ob_refcnt);
      Py_DECREF(pValueII);
      //printf("free pValueIQ\t from refcount\t %ld\n", pValueIQ->ob_refcnt);
      Py_DECREF(pValueIQ);
      //printf("free pValueQI\t from refcount\t %ld\n", pValueQI->ob_refcnt);
      Py_DECREF(pValueQI);
      //printf("free pValueQQ\t from refcount\t %ld\n", pValueQQ->ob_refcnt);
      Py_DECREF(pValueQQ);
	       
      // and the refpower
      //printf("free pValue\t from refcount\t %ld\n", pValue->ob_refcnt);
      Py_DECREF(pValue);
   
      // 4 X pList
      // These don't need to be freed since PyList_SetItem() steals the reference
      /*
      printf("free pListII\t from refcount\t %ld\n", pListII->ob_refcnt);
      Py_DECREF(pListII);
      printf("free pListIQ\t from refcount\t %ld\n", pListIQ->ob_refcnt);
      Py_DECREF(pListIQ);
      printf("free pListQI\t from refcount\t %ld\n", pListQI->ob_refcnt);
      Py_DECREF(pListQI);
      printf("free pListQQ\t from refcount\t %ld\n", pListQQ->ob_refcnt);
      Py_DECREF(pListQQ);
      */
     
      if (DEBUG)
         printf("all Python Objects from callback() freed\n\n");
 
   } // END SPECTRA LOOP
/////////////////////////////  END LOOP OVER ALL SPECTRA IN FILE  ///////////////////////////////////
   
   fclose(fp);   //close input data file
		    
   if (!error){
      //ouput stats for last spectra in file
      printf("UNIXTIME is %" PRIu64 "\n", UNIXTIME);
      printf("CORRTIME is %.6f\n", (corr.corrtime*256.)/(FS_FREQ*1000000.));
      printDateTimeFromEpoch((long long) UNIXTIME);
      printf("UNIT is %d\n", UNIT);
      printf("DEV  is %d\n",  DEV); 
      printf("NINT is %d\n", NINT);
      printf("%.2f %.2f %.2f %.2f\n",  (double)corr.Ihi/(double)corr.corrtime, \
                                       (double)corr.Qhi/(double)corr.corrtime, \
                                       (double)corr.Ilo/(double)corr.corrtime, \
                                       (double)corr.Qlo/(double)corr.corrtime);
      printf("nlags=%d\n", N);
      printf("etaQ\t%.3f\n", 1/sqrt(P_I*P_Q));
      printf("ETAQ\t%.3f\n\n", 1/sqrt(Ipwr*Qpwr));

      // current scanID, and scanID used for cal
      printf("The cal    is from scanID: %d\n", CALID);
      printf("The data   is from scanID: %d\n", scanID);
   }

   //timing
   gettimeofday(&end, 0);
   int seconds = end.tv_sec - begin.tv_sec;
   int microseconds = end.tv_usec - begin.tv_usec;
   double elapsed = seconds + microseconds*1e-6;
   if (!error)
      printf("AVG FFTW %.1f ms in %d spectra\n\n", 1000.*elapsed/(sz/bps), sz/bps);
   fflush(stdout);

   
//
// move to within 100 spectra loop
//
/*
   // Clean Up memory before leaving callback()
   // All of these objects are malloced at the start of callback but re-used every spectra
   // 4 X pArgs 
   printf("free pArgsII\t from refcount\t %ld\n", pArgsII->ob_refcnt);
   Py_DECREF(pArgsII);
   printf("free pArgsIQ\t from refcount\t %ld\n", pArgsIQ->ob_refcnt);
   Py_DECREF(pArgsIQ);
   printf("free pArgsQI\t from refcount\t %ld\n", pArgsQI->ob_refcnt);
   Py_DECREF(pArgsQI);
   printf("free pArgsQQ\t from refcount\t %ld\n", pArgsQQ->ob_refcnt);
   Py_DECREF(pArgsQQ);
    
   // and the refpower
   printf("free pArgs\t from refcount\t %ld\n", pArgs->ob_refcnt);
   Py_DECREF(pArgs);

   // 4 X pValue
   printf("free pValueII\t from refcount\t %ld\n", pValueII->ob_refcnt);
   Py_DECREF(pValueII);
   printf("free pValueIQ\t from refcount\t %ld\n", pValueIQ->ob_refcnt);
   Py_DECREF(pValueIQ);
   printf("free pValueQI\t from refcount\t %ld\n", pValueQI->ob_refcnt);
   Py_DECREF(pValueQI);
   printf("free pValueQQ\t from refcount\t %ld\n", pValueQQ->ob_refcnt);
   Py_DECREF(pValueQQ);
	       
   // and the refpower
   printf("free pValue\t from refcount\t %ld\n", pValue->ob_refcnt);
   Py_DECREF(pValue);
   

   // 4 X pList
   // These don't need to be freed since PyList_SetItem() steals the reference
   
   printf("free pListII\t from refcount\t %ld\n", pListII->ob_refcnt);
   Py_DECREF(pListII);
   printf("free pListIQ\t from refcount\t %ld\n", pListIQ->ob_refcnt);
   Py_DECREF(pListIQ);
   printf("free pListQI\t from refcount\t %ld\n", pListQI->ob_refcnt);
   Py_DECREF(pListQI);
   printf("free pListQQ\t from refcount\t %ld\n", pListQQ->ob_refcnt);
   Py_DECREF(pListQQ);
   
     
   printf("all Python Objects from callback() freed\n");
 */
 //
 // end test free in loop
 //
     
   // Free filename type char
   free(prefix);

   // Free Influx query string
   free(query);

   // Free Influx query string
   free(fullpath);
   
}






