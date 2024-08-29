#include <math.h>
#include <it/math.h>
#include <fftw3.h>
#include <sys/time.h>
#include <Python.h>
#include <stdbool.h>


// function to process correlator lag file
void const callback(char *filein);

// structure to hold lag data
struct corrType
{
   uint32_t corrtime;      //32 bit length of accumulation
   uint32_t Ihi;           //32 bit Ihi monitor
   uint32_t Qhi;           //32 bit Qhi monitor
   uint32_t Ilo;           //32 bit Ilo monitor
   uint32_t Qlo;           //32 bit Qlo monitor
   uint32_t Ierr;          //32 bit Ierr monitor
   uint32_t Qerr;          //32 bit Qerr monitor
   int32_t *II;            //32 bit II lag values
   int32_t *QI;            //32 bit QI lag values
   int32_t *IQ;            //32 bit IQ lag values
   int32_t *QQ;            //32 bit QQ lag values
   float *IIqc;            //32 bit II lag values (normlized float quantization corrected)
   float *QIqc;            //32 bit QI lag values (normlized float quantization corrected)
   float *IQqc;            //32 bit IQ lag values (normlized float quantization corrected)
   float *QQqc;            //32 bit QQ lag values (normlized float quantization corrected)
};

// structure to hold fits header pointer
typedef struct s_header {
	int		unit;
	int		dev;
	int		mixer;
	int		nint;
	double		fulltime;
	int		nbytes;
	int  		corrtime;
	float		inttime;
	uint32_t	row_flag;
	uint16_t	channel_flag;
	int		Ihi;
	int		Qhi;
	int		Ilo;
	int		Qlo;
	int		Ierr;
	int		Qerr;
	float		VIhi;
	float		VQhi;
	float		VIlo;
	float		VQlo;
	int		scanID;
	int		subScan;
	int		CALID;
	float		THOT;
	float		RA;
	float		DEC;
	float		LO;
	float		IF;
	float 		VLSR;
	char		*type;
	char		*filename;
	float           *psat;
	float           *imon;
	float           *gmon;
} s_header;


