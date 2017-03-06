/* miscellaneous scalar operations */

float abx (a)
float a;
{
	float x;

	if (a >= 0.) x =  a;
	if (a <  0.) x = -a;
	return (x);
}

float power (a, b)
float a, b;
{
	float x;

	x = b*log (a);
	x = exp (x);
	return (x);
}

float fr (a)
float a;
{
	float x;

	x = exp (-0.5*a*a);
	return (x);
}

float rms (a, X)
int a;
float X[750];
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

float norm (a, X)
int a;
float X[750];
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

float inner (a, X, Y)
int a;
float X[750], Y[750];
{
	int i;
	float x;

	x = 0.;
	for (i=1; i<=a; ++i) {
		x = x + X[i]*Y[i];
	}
	return (x);
}

void outer (a, X, Y, Z)
int a;
float X[750], Y[750], Z[750][750];
{
	int i, j;

	for (i=1; i<=a; ++i) {
		for (j=1; j<=a; ++j) {
			Z[i][j] = X[i]*Y[j];
		}
	}
}

/* matrix operations */

void AX1 (a, b, A, X, Y)
int a, b;
float A[750][750], X[750], Y[750];
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

void AX2 (a, b, A, X, Y)
int a, b;
float A[750][750], X[750], Y[750];
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

void scale (a, b, A, B)
int a;
float b, A[750][750], B[750][750];
{
	int i, j;

	for (i=1; i<=a; ++i) {
		for (j=1; j<=a; ++j) {
			B[i][j] = b*A[i][j];
		}
	}
}

void Madd (a, b, A, B, C)
int a, b;
float A[750][750], B[750][750], C[750][750];
{
	int i, j;

	for (i=1; i<=a; ++i) {
		for (j=1; j<=b; ++j) {
			C[i][j] = A[i][j] + B[i][j];
		}
	}
}

void Msub (a, b, A, B, C)
int a, b;
float A[750][750], B[750][750], C[750][750];
{
	int i, j;

	for (i=1; i<=a; ++i) {
		for (j=1; j<=b; ++j) {
			C[i][j] = A[i][j] - B[i][j];
		}
	}
}

void Mmul (a, b, c, A, B, C)
int a, b, c;
float A[750][750], B[750][750], C[750][750];
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
int a;
float s, X[10], M[10];
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

float kerny (a, s, X, M)
int a;
float s, X[10], M[10];
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
