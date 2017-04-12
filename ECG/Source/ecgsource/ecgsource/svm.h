#ifndef _LIBSVM_H
#define _LIBSVM_H

#define LIBSVM_VERSION 321

#include <iostream>
#include <fstream>
#include <cstdio>
#pragma warning(disable:4996)
using namespace std;
extern int libsvm_version;


struct svm_node
{
	int index;
	double value;
};

struct svm_problem
{
	int l;
	double *y;
	struct svm_node **x;
};

enum { C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR };	/* svm_type */
enum { LINEAR, POLY, RBF, SIGMOID, PRECOMPUTED }; /* kernel_type */

struct svm_parameter
{
	int svm_type;
	int kernel_type;
	int degree;	/* for poly */
	double gamma;	/* for poly/rbf/sigmoid */
	double coef0;	/* for poly/sigmoid */

					/* these are for training only */
	double cache_size; /* in MB */
	double eps;	/* stopping criteria */
	double C;	/* for C_SVC, EPSILON_SVR and NU_SVR */
	int nr_weight;		/* for C_SVC */
	int *weight_label;	/* for C_SVC */
	double* weight;		/* for C_SVC */
	double nu;	/* for NU_SVC, ONE_CLASS, and NU_SVR */
	double p;	/* for EPSILON_SVR */
	int shrinking;	/* use the shrinking heuristics */
	int probability; /* do probability estimates */
};

//
// svm_model
// 
struct svm_model
{
	struct svm_parameter param;	/* parameter */
	int nr_class;		/* number of classes, = 2 in regression/one class svm */
	int l;			/* total #SV */
	struct svm_node **SV;		/* SVs (SV[l]) */
	double **sv_coef;	/* coefficients for SVs in decision functions (sv_coef[k-1][l]) */
	double *rho;		/* constants in decision functions (rho[k*(k-1)/2]) */
	double *probA;		/* pariwise probability information */
	double *probB;
	int *sv_indices;        /* sv_indices[0,...,nSV-1] are values in [1,...,num_traning_data] to indicate SVs in the training set */

							/* for classification only */

	int *label;		/* label of each class (label[k]) */
	int *nSV;		/* number of SVs for each class (nSV[k]) */
					/* nSV[0] + nSV[1] + ... + nSV[k-1] = l */
					/* XXX */
	int free_sv;		/* 1 if svm_model is created by svm_load_model*/
						/* 0 if svm_model is created by svm_train */
};


struct svm_model *svm_train(const struct svm_problem *prob, const struct svm_parameter *param);
void svm_cross_validation(const struct svm_problem *prob, const struct svm_parameter *param, int nr_fold, double *target);

int svm_save_model(const char *model_file_name, const struct svm_model *model);
struct svm_model *svm_load_model(const char *model_file_name);

int svm_get_svm_type(const struct svm_model *model);
int svm_get_nr_class(const struct svm_model *model);
void svm_get_labels(const struct svm_model *model, int *label);
void svm_get_sv_indices(const struct svm_model *model, int *sv_indices);
int svm_get_nr_sv(const struct svm_model *model);
double svm_get_svr_probability(const struct svm_model *model);

double svm_predict_values(const struct svm_model *model, const struct svm_node *x, double* dec_values);
double svm_predict(const struct svm_model *model, const struct svm_node *x);
double svm_predict_probability(const struct svm_model *model, const struct svm_node *x, double* prob_estimates);

void svm_free_model_content(struct svm_model *model_ptr);
void svm_free_and_destroy_model(struct svm_model **model_ptr_ptr);
void svm_destroy_param(struct svm_parameter *param);

const char *svm_check_parameter(const struct svm_problem *prob, const struct svm_parameter *param);
int svm_check_probability_model(const struct svm_model *model);

void svm_set_print_string_function(void (*print_func)(const char *));


#ifdef __cplusplus

class SVM {
public:
	int n;
	int dim, npos, nneg;
	double **pos, **neg;
	struct svm_parameter param;
	struct svm_problem prob;
	struct svm_model *model;

	int TP, TN, FP, FN;
	double acc, pre1, pre2, rec1, rec2, f11, f12;

	void clear() {
		int i;

		for (i = 0; i < n; i++)
			free(prob.x[i]);

		free(prob.x);
		free(prob.y);
		free(model);
	}

	double getOutput(double *pat) {
		int i;
		struct svm_node *x = (struct svm_node*)malloc(sizeof(struct svm_node) * (dim + 1));

		for (i = 0; i < dim; i++) {
			x[i].index = i + 1;
			x[i].value = pat[i];
		}

		x[i].index = -1;

		double dec_values = 0.f;
		double pred = svm_predict_values(model, x, &dec_values);

		free(x);

		return dec_values;
	}

	bool isPos(double *pat) {
		int i;
		bool result;
		struct svm_node *x = (struct svm_node*)malloc(sizeof(struct svm_node) * (dim + 1));

		for (i = 0; i < dim; i++) {
			x[i].index = i + 1;
			x[i].value = pat[i];
		}

		x[i].index = -1;
		
		double dec_values = 0.f;
		double pred = svm_predict_values(model, x, &dec_values);

		result = (pred == 1 ? true : false);
		free(x);

		return result;
	}

	//Test
	void test_counting(int test_npos, int test_nneg, double **test_pos, double **test_neg) {
		int i;

		for (i = 0; i < test_npos; i++) {
			if (isPos(test_pos[i])) TP++;
			else FP++;
		}

		for (i = 0; i < test_nneg; i++) {
			if (isPos(test_neg[i])) FN++;
			else TN++;
		}
	}

	void test_measure(char *fname) {
		char fch[100];
		ofstream fout;

		sprintf(fch, "%s_svm.csv", fname);
		fout.open(fch);

		acc = 1.f - ((double)(FP + FN)) / (TP + FP + TN + FN);
		pre1 = (double)TP / (TP + FP);
		rec1 = (double)TP / (TP + FN);
		f11 = (2.f * pre1 * rec1) / (pre1 + rec1);
		pre2 = (double)TN / (TN + FN);
		rec2 = (double)TN / (TN + FP);
		f12 = (2.f * pre2 * rec2) / (pre2 + rec2);

		cout << "---------------SVM---------------" << endl;
		cout << "Accuracy : " << acc << endl;
		cout << "---------------NR----------------" << endl;
		cout << "Precision (NR) : " << pre1 << endl;
		cout << "Recall : " << rec1 << endl;
		cout << "F1 : " << f11 << endl;
		cout << "---------------AN----------------" << endl;
		cout << "Precision (NR) : " << pre2 << endl;
		cout << "Recall : " << rec2 << endl;
		cout << "F1 : " << f12 << endl;
		cout << endl;

		fout << "Accuracy," << acc << endl;
		fout << "Precision," << pre1 << endl;
		fout << "Recall," << rec1 << endl;
		fout << "F1," << f11 << endl;
		fout << "Precision," << pre2 << endl;
		fout << "Recall," << rec2 << endl;
		fout << "F1," << f12 << endl;
		fout.close();
	}

	//Training
	void train(int dim, int npos, int nneg, double **pos, double **neg) {
		int i, j;
		this->dim = dim;
		this->npos = npos;
		this->nneg = nneg;
		this->pos = pos;
		this->neg = neg;

		n = npos + nneg;
		prob.l = n;
		prob.y = (double*)malloc(sizeof(double) * n);
		prob.x = (struct svm_node**)malloc(sizeof(struct svm_node*) * n);

		for (i = 0; i < npos; i++) {
			prob.x[i] = (struct svm_node*)malloc(sizeof(struct svm_node) * (dim + 1));

			for (j = 0; j < dim; j++) {
				prob.x[i][j].index = j + 1;
				prob.x[i][j].value = pos[i][j];
			}

			prob.x[i][j].index = -1;

			prob.y[i] = 1.f;

			if (!prob.x) 
				i = i;
		}

		for (i = 0; i < nneg; i++) {
			prob.x[i + npos] = (struct svm_node*)malloc(sizeof(struct svm_node) * (dim+1));

			for (j = 0; j < dim; j++) {
				prob.x[i + npos][j].index = j + 1;
				prob.x[i + npos][j].value = neg[i][j];
			}
			
			prob.x[i+ npos][j].index = -1;
			prob.y[i + npos] = -1.f;
		}

		model = svm_train(&prob, &param);
	}


	SVM() {
		// default values
		param.svm_type = C_SVC;
		param.kernel_type = LINEAR;
		param.degree = 3;
		param.gamma = 0;	// 1/num_features
		param.coef0 = 0;
		param.nu = 0.5;
		param.cache_size = 100;
		param.C = 1;
		param.eps = 1e-3;
		param.p = 0.1;
		param.shrinking = 1;
		param.probability = 0;
		param.nr_weight = 0;
		param.weight_label = NULL;
		param.weight = NULL;
	}

	SVM(double gamma) {
		// default values
		param.svm_type = C_SVC;
		param.kernel_type = RBF;
		param.degree = 3;
		param.gamma = gamma;	// 1/num_features
		param.coef0 = 0;
		param.nu = 0.5;
		param.cache_size = 100;
		param.C = 1;
		param.eps = 1e-3;
		param.p = 0.1;
		param.shrinking = 1;
		param.probability = 0;
		param.nr_weight = 0;
		param.weight_label = NULL;
		param.weight = NULL;
	}

	SVM(int svm_type, int kernel_type) {
		// default values
		param.svm_type = svm_type;
		param.kernel_type = kernel_type;
		param.degree = 3;
		param.gamma = 0;	// 1/num_features
		param.coef0 = 0;
		param.nu = 0.5;
		param.cache_size = 100;
		param.C = 1;
		param.eps = 1e-3;
		param.p = 0.1;
		param.shrinking = 1;
		param.probability = 0;
		param.nr_weight = 0;
		param.weight_label = NULL;
		param.weight = NULL;
	}

	SVM(struct svm_parameter p) {
		param.svm_type = p.svm_type;
		param.kernel_type = p.kernel_type;
		param.degree = p.degree;
		param.gamma = p.gamma;	// 1/num_features
		param.coef0 = p.coef0;
		param.nu = p.nu;
		param.cache_size = p.cache_size;
		param.C = p.C;
		param.eps = p.eps;
		param.p = p.p;
		param.shrinking = p.shrinking;
		param.probability = p.probability;
		param.nr_weight = p.nr_weight;
		param.weight_label = p.weight_label;
		param.weight = p.weight;
	}

};

#endif


#endif /* _LIBSVM_H */
