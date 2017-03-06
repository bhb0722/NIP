#define di1 800
#define id1 10
/* miscellaneous scalar operations */
/* absoulte value */
float abx (a)
float a;
{
	float x;

	if (a >= 0.) x =  a;
	if (a <  0.) x = -a;
	return (x);
}
/* a^b */
float power (a, b)
float a, b;
{
	float x;

	x = b*log (a);
	x = exp (x);
	return (x);
}
/* value of the function exp(-(1/2)*x^2 */ 
float fr (a)
float a;
{
	float x;

	x = exp (-0.5*a*a);
	return (x);
}
/* root mean square */
float rms (a, X)
int a; /* # of size */
float X[di1]; /* array */
{
	int i;
	float x;

	x = 0;
	for (i=1; i<=a; ++i) {
		x = x + X[i]*X[i];
	}
	x = x/a;
	x = sqrt (x);
	return (x);
}

/* vector operations */
/* Euclidean Norm */
float norm (a, X)
int a; /* # of size */
float X[di1]; /* array */
{
	int i;
	float x;

	x = 0.;
	for (i=1; i<=a; ++i) {
		x = x + X[i]*X[i];
	}
	x = sqrt (x);
	return (x);
}
/* X.Y */
float inner (a, X, Y)
int a; /* # of size */
float X[di1], Y[di1]; /* array */
{
	int i;
	float x;

	x = 0.;
	for (i=1; i<=a; ++i) {
		x = x + X[i]*Y[i];
	}
	return (x);
}
/* c[i][j] = a[i]*b[j] */
void outer (a, X, Y, Z)
int a; /* # of size */
float X[di1], Y[di1], Z[di1][di1]; /* inputs, output */
{
	int i, j;

	for (i=1; i<=a; ++i) {
		for (j=1; j<=a; ++j) {
			Z[i][j] = X[i]*Y[j];
		}
	}
}

/* matrix operations */
/* Y = AX; */
void AX1 (a, b, A, X, Y)
int a, b; /* the size of rows, the size of columns */
float A[di1][di1], X[di1], Y[di1]; /* inputs , output */
{
	int i, j;
	float x;

	for (i=1; i<=a; ++i) {
		Y[i] = 0.;
		for (j=1; j<=b; ++j) {
			Y[i] = Y[i] + A[i][j]*X[j];
		}
	}
}
/* Y = A'X */
void AX2 (a, b, A, X, Y)
int a, b; /* the size of rows, the size of columns */
float A[di1][di1], X[di1], Y[di1]; /* inputs, output */
{
	int i, j;
	float x;

	for (i=1; i<=b; ++i) {
		Y[i] = 0.;
		for (j=1; j<=a; ++j) {
			Y[i] = Y[i] + A[j][i]*X[j];
		}
	}
}
/* B = aA */
void scale (a, b, A, B)
int a; /* the size of matrix A */
float b, A[di1][di1], B[di1][di1]; /* b : scale factor, A = input, B = output */
{
	int i, j;

	for (i=1; i<=a; ++i) {
		for (j=1; j<=a; ++j) {
			B[i][j] = b*A[i][j];
		}
	}
}
/* Matrix add : C = A + B */
void Madd (a, b, A, B, C)
int a, b; /* the size of row, the size of column of matrix */
float A[di1][di1], B[di1][di1], C[di1][di1]; /* matrices */
{
	int i, j;

	for (i=1; i<=a; ++i) {
		for (j=1; j<=b; ++j) {
			C[i][j] = A[i][j] + B[i][j];
		}
	}
}
/* Matrix sub : C = A - B */
void Msub (a, b, A, B, C)
int a, b;
float A[di1][di1], B[di1][di1], C[di1][di1];
{
	int i, j;

	for (i=1; i<=a; ++i) {
		for (j=1; j<=b; ++j) {
			C[i][j] = A[i][j] - B[i][j];
		}
	}
}
/* C = A*B */
void Mmul (a, b, c, A, B, C)
int a, b, c;
float A[di1][di1], B[di1][di1], C[di1][di1];
{
	int i, j, k;
	float x;

	for (i=1; i<=a; ++i) {
		for (j=1; j<=c; ++j) {
			x = 0.;
			for (k=1; k<=b; ++k) {
				x = x + A[i][k]*B[k][j];
			}
			C[i][j] = x;
		}
	}
}

/* chaotic time-series */

float LM (a)
float a;
{
	float x;

	x = 4*a*(1-a);
	return (x);
}

float MG (a, b)
float a, b;
{
	float x, y, z;

	x = 10*log (a);
	x = exp (x);
	y = 0.2*a/(1 + x);
	z = y - 0.1*b;
	return (z);
}

/* Gaussian kernel functions */

float kernx (a, s, X, M)
int a; /* the size of X, M */
float s, X[id1], M[id1]; /* s = sigma ; X = input ; M = mean */
{
	int i;
	float x, y, z;

	z = 1;
	x = 0;
	for (i=1; i<=a; ++i) {
		x = x + (X[i]-M[i])*(X[i]-M[i]);
	}
	y = -x/(2*s*s);
	z = z*exp (y);
	return (z);
}
/* the derivative value of Gaussian Kernel Function */
float kerny (a, s, X, M)
int a; /* size of vector */
float s, X[id1], M[id1]; /* s = sigma ; X = input ; M = mean */
{
	int i;
	float x, y, z;

	z = 1;
	x = 0;
	for (i=1; i<=a; ++i) {
		x = x + (X[i]-M[i])*(X[i]-M[i]);
	}
	y = -x*s;
	z = z*exp (y);
	return (z);
}
