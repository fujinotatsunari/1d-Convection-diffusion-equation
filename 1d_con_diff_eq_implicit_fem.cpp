/****************************************************************************/
//   Galerkin�L���v�f�@�ɂ��ꎟ���ڗ��g�U�������̐��l�v�Z�V�~�����[�V����(���S�A��@)
/****************************************************************************/
#include <cstdio>
#include<stdlib.h>
#include <cstdio>
#include <cmath>
#define _USE_MATH_DEFINES
#include<math.h>
#pragma warning(disable:4996)

#define NODE 2


void gaussian_elimination(int dim, double** mat, double* b, double* x);


int main(void)
{
	double x_left = 0.0, x_right = 1.0; // �v�Z�̈�̎w��
	int nelement = 40; // �v�f��
	int dim = nelement + 1;// �s��A�̗v�f��
	int nnode = nelement + 1; // �ߓ_��
	double alpha = 0.1; // �g�U�W��
	double u = 0.5;//����
	double dt = 1.0e-2; // ���ԍ��ݕ�
	int nend = 1000; // ���ԃX�e�b�v����

	/*********** �ߓ_���W�̌v�Z ***********/
	double dx = (x_right - x_left) / (double)nelement; // �i�q��
	double* x = new double[nnode];

	for (int i = 0; i < nnode; i++) {
		x[i] = x_left + dx * (double)i; // x���W
	}



	/*********** ���[�J���Ȑߓ_�ԍ��ƃO���[�o���Ȑߓ_�ԍ��̑Ή��t�� ***********/
	int** nbool;
	nbool = new int* [nelement];
	nbool[0] = new int[nelement * NODE];
	for (int ie = 0; ie < nelement; ie++) {
		nbool[ie] = nbool[0] + ie * NODE;
	}

	for (int ie = 0; ie < nelement; ie++) {
		nbool[ie][0] = ie; // �v�f�ԍ�ie�̋Ǐ��ߓ_�ԍ�0�ɑΉ�����O���[�o���Ȑߓ_�ԍ�
		nbool[ie][1] = ie + 1;
	}

	/*********** ���R�r�A���̌v�Z ***********/
	double* jacobian = new double[nelement];

	for (int ie = 0; ie < nelement; ie++) {
		int i1 = nbool[ie][0];
		int i2 = nbool[ie][1];

		jacobian[ie] = 0.5 * (x[i2] - x[i1]);
	}
	//�v�Z��ԏ�̋Ǐ��ߓ_�ԍ�
	double xi[NODE] = { -1.0, 1.0 };

	/*********** �W�����ʍs�� ***********/
	double* fm = new double[nnode];

	for (int i = 0; i < nnode; i++) {
		fm[i] = 0.0;
	}

	for (int j = 0; j < NODE; j++) {
		for (int ie = 0; ie < nelement; ie++) {
			int nu = nbool[ie][j];
			fm[nu] += jacobian[ie];
		}
	}

	/********** ���ʍs�� *********/

	double em[NODE][NODE];

	for (int i = 0; i < NODE; i++) {
		for (int j = 0; j < NODE; j++) {
			em[i][j] = 0.5 * (1 + (xi[i] * xi[j]) / 3.0);//���R�r�A���Ȃ�

		}
	}

	/*********** �v�f�g�U�s�� ***********/

	double dm[NODE][NODE];

	for (int i = 0; i < NODE; i++) {
		for (int j = 0; j < NODE; j++) {
			dm[i][j] = 0.5 * xi[i] * xi[j];//���R�r�A���Ȃ�

		}

	}

	/*********** �v�f�ڗ��s�� ***********/

	double am[NODE][NODE];

	for (int i = 0; i < NODE; i++) {
		for (int j = 0; j < NODE; j++) {
			am[i][j] = 0.5 * xi[j];//���R�r�A���Ȃ�

		}

	}



	/*********  �W���s��̘a(���Ӎs��̗v�f) (���[Dirichlet) **********/
	double* sm1 = new double[nnode];
	double* sm2 = new double[nnode];
	double* sm3 = new double[nnode];


	for (int i = 0; i < nnode; i++) {
		if (i == 0) {//���[
			sm1[i] = 0;
			sm2[i] = 0;
			sm3[i] = 0;

		}
		else if (i == nnode-1) {//�E�[
			sm1[i] = 0;
			sm2[i] = 0;
			sm3[i] = 0;
		}
		else {
			sm1[i] = em[1][0] * jacobian[i - 1] + u * dt * am[1][0] + alpha * dt * dm[1][0] / jacobian[i - 1];
			sm2[i] = em[0][0] * jacobian[i - 1] + u * dt * am[0][0] + alpha * dt * dm[0][0] / jacobian[i - 1] + em[1][1] * jacobian[i] + u * dt * am[1][1] + alpha * dt * dm[1][1] / jacobian[i];
			sm3[i] = em[0][1] * jacobian[i] + u * dt * am[0][1] + alpha * dt * dm[0][1] / jacobian[i];

		}

	}
	/*for (int i = 0; i < nnode; i++) {
		printf("sm1=%23.15e: sm2=%23.15e: sm3=%23.15e\n", sm1[i], sm2[i], sm3[i]);
	}*/

	/*********** ���������̐ݒ� ***********/
	double* theta = new double[nnode];

	for (int i = 0; i < nnode; i++) {
		if (x[i] >= 0.4 && x[i] <= 0.6 ) {
			theta[i] = 1.0;
		}
		else {
			theta[i] = 0.0;
		}
	}




	/************ �s��������̐ݒ�  ********/
	double** Amat; // ���Ӎs�� A
	double* phi_b; // �E�Ӄx�N�g��b ���݂̃X�e�b�v�̃ӂ���Ȃ�
	double* phi_x; // ���x�N�g��x ���X�e�b�v�̃�
	double** phi_s; //��ۑ��p


	/* �������̓��I�m�� */
	Amat = new double* [dim];
	phi_b = new double [dim];
	phi_x = new double [dim];
	phi_s = new double* [dim];

	for (int i = 0; i < dim; i++) {
		Amat[i] = new double[dim];
	}
	for (int i = 0; i < nnode; i++) {
		phi_s[i] = new double[nend];
	}

	/* ������ */
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			Amat[i][j] = 0.0;
		}
		phi_b[i] = 0.0;
		phi_x[i] = 0.0;

	}

	for (int i = 0; i < nnode; i++) {
		for (int j = 0; j < nend; j++) {
			phi_s[i][j] = 0.0;
		}
	}



	/*�s��̗v�f����(���[DirichletB.C)*/
	for (int i = 0; i < nnode; i++) {
		for (int j = 0; j < nnode; j++) {
			if (i == 0) {
				Amat[i][i] = 1.0;

			}
			else if (i == nnode - 1) {
				Amat[i][i] = 1.0;
			}
			else {
				Amat[i][i - 1] = sm1[i];
				Amat[i][i] = sm2[i];
				Amat[i][i + 1] = sm3[i];
			}
		}
	}

	
	/*���Ӄx�N�g��b�̗v�f���*/
	double dirichlet_left = 0.0;//dirichlet���E�����l
	double dirichlet_right = 0.0;
	double neumann_h_left = 0.0;//neumann���E�����l:: h = ��*�݃�/��x
	double neumann_h_right = 0.0;
	double neumann_left = neumann_h_left * dx / alpha;
	double neumann_right = neumann_h_right * dx / alpha;

	for (int i = 0; i < nnode; i++) {
		if (i == 0) {
			phi_b[i] = dirichlet_left;
		}
		else if (i == nnode - 1) {
			phi_b[i] = dirichlet_right;
		}
		else {
			phi_b[i] = em[1][0] * jacobian[i - 1] * theta[i - 1] + (em[0][0] * jacobian[i - 1] + em[1][1] * jacobian[i]) * theta[i] + em[0][1] * jacobian[i] * theta[i + 1];
		}
	}

	/*
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			printf(" A[%d][%d] = %9f", i, j, Amat[i][j]);
		}
		printf(" b[%d] = %9f\n", i, phi_b[i]);
	}
	printf("\n");

	*/


	/**************   ���l�v�Z�̊J�n    ************/

	double time = 0.0;


	/*�������̕ۑ�*/
	for (int i = 0; i < nnode; i++) {
		phi_s[i][0] = theta[i];
	}

	
	for (int n = 1; n <= nend; n++) {
		time = (double)n * dt;//���݂̎���


		/*�������犮�S�A��@�ɂ���ĉ���*/

		/*�s�������������*/

		gaussian_elimination(dim, Amat, phi_b, phi_x);//�K�E�X�̏����@



		for (int i = 0; i < nnode; i++) {//���x�N�g���̕ۑ�

			theta[i] = phi_x[i];
			phi_s[i][n] = phi_x[i];

		}
		for (int i = 0; i < nnode; i++) {
			if (i == 0) {
				phi_b[i]= dirichlet_left;
			}
			else if (i == nnode - 1) {
				phi_b[i] = dirichlet_right;

			}
			else {
				phi_b[i] = em[1][0] * jacobian[i - 1] * theta[i - 1] + (em[0][0] * jacobian[i - 1] + em[1][1] * jacobian[i]) * theta[i] + em[0][1] * jacobian[i] * theta[i + 1];
			}
		}


		/*�s����č\��(���[DirichletB.C)*/
		for (int i = 0; i < nnode; i++) {
			for (int j = 0; j < nnode; j++) {
				if (i == 0) {
					Amat[i][i] = 1.0;

				}
				else if (i == nnode - 1) {
					Amat[i][i] = 1.0;
				}
				else {
					Amat[i][i - 1] = sm1[i];
					Amat[i][i] = sm2[i];
					Amat[i][i + 1] = sm3[i];
				}
			}
		}

	

	}

	
	/*********** ��ʂւ̏o�� ***********/
	for (int i = 0; i < nnode; i++) {
		printf("%23.15e %23.15e\n", x[i], theta[i]);
	}

	double *T, *X;//���Ԃƍ��W�̋L�^
	T = new double[nend];
	X = new double[nelement];
	FILE* fp;

	//csv�o��

	fp = fopen("result.csv", "w");
	fprintf(fp, ",");
	for (int n = 0; n < nend; n++) {
		T[n] = dt * n;
		fprintf(fp, "%f,", T[n]);
	}
	fprintf(fp, "\n");

	for (int i = 0; i < nnode; i++) {
		X[i] = dx * i;
		fprintf(fp, "%f,", X[i]);
		for (int n = 0; n < nend; n++) {
			fprintf(fp, "%23.15e,", phi_s[i][n]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	//���������
	for (int i = 0; i < dim; i++) {
		delete[] Amat[i];
	}
	for (int i = 0; i < nnode; i++) {
		delete[] phi_s[i];
	}
	delete[] phi_s;
	delete[] Amat;
	delete[] phi_b;
	delete[] phi_x;

	delete[] x;
	delete[] nbool[0]; delete[] nbool;
	delete[] jacobian;
	delete[] fm;
	delete[] theta;

	return 0;
}

void gaussian_elimination(int dim, double** mat, double* b, double* x) {
	/****************************************/
	/************ �K�E�X�̏����@ ************/
	/****************************************/
	// �O�i����
	for (int i = 0; i < dim - 1; i++) {
		/******************************/
		/******** �s�{�b�g�I�� ********/
		/******************************/
		int id = i;
		double gg = fabs(mat[i][i]);
		// �ő�l�̍s��T��
		for (int ip = i + 1; ip < dim; ip++) {
			double ah = fabs(mat[ip][i]);
			if (ah > gg) {
				id = ip;
				gg = ah;
			}
		}
		// �[���������
		if (gg < 1.0e-8) {
			printf("�Ίp����������������\n");
			exit(1); // �����I��
		}
		// �s�̓���ւ�
		for (int jp = i; jp < dim; jp++) {
			double w = mat[i][jp];
			mat[i][jp] = mat[id][jp];
			mat[id][jp] = w;
		}
		double w = b[i];
		b[i] = b[id];
		b[id] = w;

		/******************************/
		/* �Ίp������艺���[���ɂ��� */
		/******************************/
		double ppp = mat[i][i];
		for (int ii = i + 1; ii < dim; ii++) {
			double q = mat[ii][i] / ppp;
			for (int j = i; j < dim; j++) {
				mat[ii][j] = mat[ii][j] - q * mat[i][j];
			}
			b[ii] = b[ii] - q * b[i];
		}

		// �s��ƉE�Ӄx�N�g����\��
		/*printf("***** %d�s�� *****\n", i);
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				printf(" A[%d][%d] = %9f", i, j, mat[i][j]);
			}
			printf(" b[%d] = %9f\n", i, b[i]);
		}
		printf("\n");*/
	}

	// ��ޑ��
	x[dim - 1] = b[dim - 1] / mat[dim - 1][dim - 1];
	for (int i = dim - 2; i >= 0; i--) {
		double s = b[i];
		for (int j = i + 1; j < dim; j++) {
			s = s - mat[i][j] * x[j];
		}
		x[i] = s / mat[i][i];
	}
	
	// ���x�N�g����\��
	/*printf("***** �� *****\n");
    for(int i = 0; i < dim; i++) {
        printf(" x[ %d ] = %23.15e\n", i, x[i]);
    }*/
}