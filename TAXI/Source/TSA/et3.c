/* Identify Time Series: */
/* search for the proper time delay and embedding dimension */

#include <stdio.h>
#include <math.h>
#include "utility.h"

#define pi 3.141596
#define id 20       /* max. dimension of input space */   
#define di 750      /* max. dimension of GKFs */
#define Nm 10000    /* max. number of time series data */      
#define TkL 2       /* lower bound of embedding dimension */
#define TkU 10      /* upper bound of embedding dimension */
#define TdL 1       /* lower bound of delay time */
#define TdU 20      /* upper bound of delay time */
#define DTd 1       /* delay time step */

main(argc, argv)
int argc; char *argv[];

{

	FILE  *ifp, *ofp;
        char  ch;
	int   i, j, k, l, n, m, p, q;
	int   M, N, Nt, Nf, N1, N2, Np, Ns, Td, Tk, Tp, St, ind[30];
	float u, v, w, x, y, z;
	float ds[Nm], ts[Nm], tsi[Nm][id], tso[Nm], tmp[30];
        float etm[30][30], ets[30][30], cri[30][30];
        float e50[30][30], e75[30][30], e90[30][30];

	/* ---------- READ INPUTS ---------- */

	/* prepare test data */

        printf ("prediction time = ");
        scanf ("%d", &M);

	ifp = fopen (argv[1], "r");
	for (i=1; (ch=getc(ifp)) != EOF; ++i) {
		fscanf (ifp, "%f", &ts[i]);
	}
	close (ifp);
	N = i - 1;
        printf ("%d %d\n", M, N);
        for (Tk=TkL; Tk<=TkU; ++Tk)
        for (Td=TdL; Td<=TdU; Td=Td+DTd) {

	Ns = N - Td*(Tk-1);
	for (i=1; i<=Tk; ++i) {
		k = Td*(i-1);
		for (j=1; j<=Ns; ++j) {
			k = k + 1;
			tsi[j][i] = ts[k];
		}
	}
        Tp = M; /* prediction time */
	Ns = Ns - Tp;
	St = Td*(Tk-1) + Tp;
	k  = St;
	for (i=1; i<Ns; ++i) {
		k = k + 1;
		tso[i] = ts[k];
	}

        /* calculation of the smoothness of mapping */
        /* from the relationship between points of minimum distance */
        /* mean */
        z = 0;
        for (i=1; i<Ns; ++i) {
                u = 10000.;
                for (j=1; j<Ns; ++j)
                        if (i != j) {
                                y = 0;
                                for (k=1; k<=Tk; ++k) {
                                        x = tsi[i][k]-tsi[j][k];
                                        y = y + x*x;
                                }
                                if (y < u) {
                                        u = y;
                                        l = j;
				}
 		        }
                x = sqrt(u); 
                ds[i] = abx(tso[i]-tso[l])/x;
                z = z + ds[i];
        }
        Ns = Ns - 1;
        x = z/Ns;
        etm[Tk][Td] = x;
/*
        for(i=1; i<=Ns; i++) printf("%d:%f ", i, ds[i]);
        printf("\n");
*/
        for (i=0; i<Ns-1; i++)
        for (j=2; j<Ns-i+1; j++) 
             if (ds[j-1]>ds[j]) {
                  u = ds[j];
                  ds[j] = ds[j-1];
                  ds[j-1] = u;
	     }
/*
        for(i=1; i<=Ns; i++) printf("%d:%f ", i, ds[i]);
        printf("\n");
*/
        if (Ns%2==0) {
             e50[Tk][Td] = (ds[(int)(Ns*0.5)]+ds[(int)(Ns*0.5)+1])/2;
             e75[Tk][Td] = (ds[(int)(Ns*0.75)]+ds[(int)(Ns*0.75)+1])/2;
             e90[Tk][Td] = (ds[(int)(Ns*0.9)]+ds[(int)(Ns*0.9)+1])/2;
	     }
        else {
             e50[Tk][Td] = ds[(int)(Ns*0.5)+1];
             e75[Tk][Td] = ds[(int)(Ns*0.75)+1];
             e90[Tk][Td] = ds[(int)(Ns*0.9)+1];
             }             

        /* standard deviation */
        y = 0;
        for (i=1; i<=Ns; ++i) {
                z = ds[i] - x;
                y = y + z*z;
        }
        y = sqrt(y/(Ns-1));
        ets[Tk][Td] = y; 

        printf ("E = %2d, T = %2d\n", Tk, Td);
        }

        ofp = fopen (argv[2], "w");
        fprintf(ofp, "prediction step = %d\n", M);
        fprintf(ofp, "\n");
        fprintf(ofp, "*** mean ***\n");
        fprintf(ofp, "\n");

        for (Tk=TkL; Tk<=TkU; ++Tk) {
        fprintf (ofp, "%2d ", Tk);
        for (Td=TdL; Td<=TdU; Td=Td+DTd)
                fprintf (ofp, "%6.4f ", etm[Tk][Td]);
        fprintf (ofp, "\n");
        }

        fprintf(ofp, "\n");
        fprintf(ofp, "*** sorted mean ***\n");
        fprintf(ofp, "\n");
        for (Tk=TkL; Tk<=TkU; ++Tk) {
        fprintf (ofp, "%2d ", Tk);
        for (Td=TdL; Td<=TdU; Td=Td+DTd) {
                tmp[Td] = etm[Tk][Td];
                ind[Td] = Td;
        }
        for (i=TdL; i<=TdU; i+=DTd)
        for (j=TdL+DTd; j<=TdU; j+=DTd)
             if (tmp[j-1]>tmp[j]) {
                  x = tmp[j];
                  tmp[j] = tmp[j-1];
                  tmp[j-1] = x;
                  k = ind[j];
                  ind[j] = ind[j-1];
                  ind[j-1] = k;
	     }
        for (Td=TdL; Td<=TdU; Td=Td+DTd)
                  fprintf (ofp, "%d:%6.4f ", ind[Td], tmp[Td]);
        fprintf (ofp, "\n");
        }

        fprintf(ofp, "\n");
        fprintf(ofp, "*** standard deviation ***\n");
        fprintf(ofp, "\n");
        for (Tk=TkL; Tk<=TkU; ++Tk) {
        fprintf (ofp, "%2d ", Tk);
        for (Td=TdL; Td<=TdU; Td=Td+DTd)
                fprintf (ofp, "%6.4f ", ets[Tk][Td]);
        fprintf (ofp, "\n");
        }

        fprintf(ofp, "\n");
        fprintf(ofp, "*** sorted standard deviation ***\n");
        fprintf(ofp, "\n");
        for (Tk=TkL; Tk<=TkU; ++Tk) {
        fprintf (ofp, "%2d ", Tk);
        for (Td=TdL; Td<=TdU; Td=Td+DTd) {
                tmp[Td] = ets[Tk][Td];
                ind[Td] = Td;
        }
        for (i=TdL; i<=TdU; i+=DTd)
        for (j=TdL+DTd; j<=TdU; j+=DTd)
             if (tmp[j-1]>tmp[j]) {
                  x = tmp[j];
                  tmp[j] = tmp[j-1];
                  tmp[j-1] = x;
                  k = ind[j];
                  ind[j] = ind[j-1];
                  ind[j-1] = k;
	     }
        for (Td=TdL; Td<=TdU; Td=Td+DTd)
                  fprintf (ofp, "%d:%6.4f ", ind[Td], tmp[Td]);
        fprintf (ofp, "\n");
        }

        fprintf(ofp, "\n");
        fprintf(ofp, "*** test criteria ***\n");
        fprintf(ofp, "\n");
        for (Tk=TkL; Tk<=TkU; ++Tk) {
        fprintf (ofp, "%2d ", Tk);
        for (Td=TdL; Td<=TdU; Td=Td+DTd)
                fprintf (ofp, "%6.4f ", cri[Tk][Td]=e90[Tk][Td]);
        fprintf (ofp, "\n");
        }

        fprintf(ofp, "\n");
        fprintf(ofp, "*** sorted test criteria ***\n");
        fprintf(ofp, "\n");
        for (Tk=TkL; Tk<=TkU; ++Tk) {
        fprintf (ofp, "%2d ", Tk);
        for (Td=TdL; Td<=TdU; Td=Td+DTd) {
                tmp[Td] = cri[Tk][Td];
                ind[Td] = Td;
        }
        for (i=TdL; i<=TdU; i+=DTd)
        for (j=TdL+DTd; j<=TdU; j+=DTd)
             if (tmp[j-1]>tmp[j]) {
                  x = tmp[j];
                  tmp[j] = tmp[j-1];
                  tmp[j-1] = x;
                  k = ind[j];
                  ind[j] = ind[j-1];
                  ind[j-1] = k;
	     }
        for (Td=TdL; Td<=TdU; Td=Td+DTd)
                  fprintf (ofp, "%d:%6.4f ", ind[Td], tmp[Td]);
        fprintf (ofp, "\n");
        }
        close (ofp);
      
        ofp = fopen(argv[3], "w");
        for (Tk=1; Tk<=TkU; ++Tk) {
             if (Tk==1) {
                  fprintf (ofp, "{{");
                  for (Td=TdL; Td<=TdU; Td=Td+DTd)
                       fprintf (ofp, "%6.4f,", 0);
                  fseek (ofp, -1, 1);
                  fprintf (ofp, "},");
	     }
             else {
                  fprintf (ofp, "{");            
                  for (Td=TdL; Td<=TdU; Td=Td+DTd)
                       fprintf (ofp, "%6.4f,", etm[Tk][Td]);
                  fseek (ofp, -1, 1);
                  fprintf (ofp, "},");
	     }
        }
        fseek (ofp, -1, 1);
        fprintf (ofp, "}");
        close (ofp);     
      
        ofp = fopen(argv[4], "w");
        for (Tk=1; Tk<=TkU; ++Tk) {
             if (Tk==1) {
                  fprintf (ofp, "{{");
                  for (Td=TdL; Td<=TdU; Td=Td+DTd)
                       fprintf (ofp, "%6.4f,", 0);
                  fseek (ofp, -1, 1);
                  fprintf (ofp, "},");
	     }
             else {
                  fprintf (ofp, "{");            
                  for (Td=TdL; Td<=TdU; Td=Td+DTd)
                       fprintf (ofp, "%6.4f,", ets[Tk][Td]);
                  fseek (ofp, -1, 1);
                  fprintf (ofp, "},");
	     }
        }
        fseek (ofp, -1, 1);
        fprintf (ofp, "}");
        close (ofp);     
      
        ofp = fopen(argv[5], "w");
        for (Tk=1; Tk<=TkU; ++Tk) {
             if (Tk==1) {
                  fprintf (ofp, "{{");
                  for (Td=TdL; Td<=TdU; Td=Td+DTd)
                       fprintf (ofp, "%6.4f,", 0);
                  fseek (ofp, -1, 1);
                  fprintf (ofp, "},");
	     }
             else {
                  fprintf (ofp, "{");            
                  for (Td=TdL; Td<=TdU; Td=Td+DTd)
                       fprintf (ofp, "%6.4f,", cri[Tk][Td]);
                  fseek (ofp, -1, 1);
                  fprintf (ofp, "},");
	     }
        }
        fseek (ofp, -1, 1);
        fprintf (ofp, "}");
        close (ofp);     

}
