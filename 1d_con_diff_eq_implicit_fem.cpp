/****************************************************************************/
//   Galerkin有限要素法による一次元移流拡散方程式の数値計算シミュレーション(完全陰解法)
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
	double x_left = 0.0, x_right = 1.0; // 計算領域の指定
	int nelement = 40; // 要素数
	int dim = nelement + 1;// 行列Aの要素数
	int nnode = nelement + 1; // 節点数
	double alpha = 0.1; // 拡散係数
	double u = 0.5;//流速
	double dt = 1.0e-2; // 時間刻み幅
	int nend = 1000; // 時間ステップ総数

	/*********** 節点座標の計算 ***********/
	double dx = (x_right - x_left) / (double)nelement; // 格子幅
	double* x = new double[nnode];

	for (int i = 0; i < nnode; i++) {
		x[i] = x_left + dx * (double)i; // x座標
	}



	/*********** ローカルな節点番号とグローバルな節点番号の対応付け ***********/
	int** nbool;
	nbool = new int* [nelement];
	nbool[0] = new int[nelement * NODE];
	for (int ie = 0; ie < nelement; ie++) {
		nbool[ie] = nbool[0] + ie * NODE;
	}

	for (int ie = 0; ie < nelement; ie++) {
		nbool[ie][0] = ie; // 要素番号ieの局所節点番号0に対応するグローバルな節点番号
		nbool[ie][1] = ie + 1;
	}

	/*********** ヤコビアンの計算 ***********/
	double* jacobian = new double[nelement];

	for (int ie = 0; ie < nelement; ie++) {
		int i1 = nbool[ie][0];
		int i2 = nbool[ie][1];

		jacobian[ie] = 0.5 * (x[i2] - x[i1]);
	}
	//計算空間上の局所節点番号
	double xi[NODE] = { -1.0, 1.0 };

	/*********** 集中質量行列 ***********/
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

	/********** 質量行列 *********/

	double em[NODE][NODE];

	for (int i = 0; i < NODE; i++) {
		for (int j = 0; j < NODE; j++) {
			em[i][j] = 0.5 * (1 + (xi[i] * xi[j]) / 3.0);//ヤコビアンなし

		}
	}

	/*********** 要素拡散行列 ***********/

	double dm[NODE][NODE];

	for (int i = 0; i < NODE; i++) {
		for (int j = 0; j < NODE; j++) {
			dm[i][j] = 0.5 * xi[i] * xi[j];//ヤコビアンなし

		}

	}

	/*********** 要素移流行列 ***********/

	double am[NODE][NODE];

	for (int i = 0; i < NODE; i++) {
		for (int j = 0; j < NODE; j++) {
			am[i][j] = 0.5 * xi[j];//ヤコビアンなし

		}

	}



	/*********  係数行列の和(左辺行列の要素) (両端Dirichlet) **********/
	double* sm1 = new double[nnode];
	double* sm2 = new double[nnode];
	double* sm3 = new double[nnode];


	for (int i = 0; i < nnode; i++) {
		if (i == 0) {//左端
			sm1[i] = 0;
			sm2[i] = 0;
			sm3[i] = 0;

		}
		else if (i == nnode-1) {//右端
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

	/*********** 初期条件の設定 ***********/
	double* theta = new double[nnode];

	for (int i = 0; i < nnode; i++) {
		if (x[i] >= 0.4 && x[i] <= 0.6 ) {
			theta[i] = 1.0;
		}
		else {
			theta[i] = 0.0;
		}
	}




	/************ 行列方程式の設定  ********/
	double** Amat; // 左辺行列 A
	double* phi_b; // 右辺ベクトルb 現在のステップのφからなる
	double* phi_x; // 解ベクトルx 次ステップのφ
	double** phi_s; //解保存用


	/* メモリの動的確保 */
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

	/* 初期化 */
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



	/*行列の要素を代入(両端DirichletB.C)*/
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

	
	/*左辺ベクトルbの要素代入*/
	double dirichlet_left = 0.0;//dirichlet境界条件値
	double dirichlet_right = 0.0;
	double neumann_h_left = 0.0;//neumann境界条件値:: h = α*∂φ/∂x
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


	/**************   数値計算の開始    ************/

	double time = 0.0;


	/*初期解の保存*/
	for (int i = 0; i < nnode; i++) {
		phi_s[i][0] = theta[i];
	}

	
	for (int n = 1; n <= nend; n++) {
		time = (double)n * dt;//現在の時刻


		/*ここから完全陰解法によって解く*/

		/*行列方程式を解く*/

		gaussian_elimination(dim, Amat, phi_b, phi_x);//ガウスの消去法



		for (int i = 0; i < nnode; i++) {//解ベクトルの保存

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


		/*行列を再構成(両端DirichletB.C)*/
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

	
	/*********** 画面への出力 ***********/
	for (int i = 0; i < nnode; i++) {
		printf("%23.15e %23.15e\n", x[i], theta[i]);
	}

	double *T, *X;//時間と座標の記録
	T = new double[nend];
	X = new double[nelement];
	FILE* fp;

	//csv出力

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

	//メモリ解放
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
	/************ ガウスの消去法 ************/
	/****************************************/
	// 前進消去
	for (int i = 0; i < dim - 1; i++) {
		/******************************/
		/******** ピボット選択 ********/
		/******************************/
		int id = i;
		double gg = fabs(mat[i][i]);
		// 最大値の行を探索
		for (int ip = i + 1; ip < dim; ip++) {
			double ah = fabs(mat[ip][i]);
			if (ah > gg) {
				id = ip;
				gg = ah;
			}
		}
		// ゼロ割を回避
		if (gg < 1.0e-8) {
			printf("対角成分が小さすぎる\n");
			exit(1); // 強制終了
		}
		// 行の入れ替え
		for (int jp = i; jp < dim; jp++) {
			double w = mat[i][jp];
			mat[i][jp] = mat[id][jp];
			mat[id][jp] = w;
		}
		double w = b[i];
		b[i] = b[id];
		b[id] = w;

		/******************************/
		/* 対角成分より下をゼロにする */
		/******************************/
		double ppp = mat[i][i];
		for (int ii = i + 1; ii < dim; ii++) {
			double q = mat[ii][i] / ppp;
			for (int j = i; j < dim; j++) {
				mat[ii][j] = mat[ii][j] - q * mat[i][j];
			}
			b[ii] = b[ii] - q * b[i];
		}

		// 行列と右辺ベクトルを表示
		/*printf("***** %d行目 *****\n", i);
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				printf(" A[%d][%d] = %9f", i, j, mat[i][j]);
			}
			printf(" b[%d] = %9f\n", i, b[i]);
		}
		printf("\n");*/
	}

	// 後退代入
	x[dim - 1] = b[dim - 1] / mat[dim - 1][dim - 1];
	for (int i = dim - 2; i >= 0; i--) {
		double s = b[i];
		for (int j = i + 1; j < dim; j++) {
			s = s - mat[i][j] * x[j];
		}
		x[i] = s / mat[i][i];
	}
	
	// 解ベクトルを表示
	/*printf("***** 解 *****\n");
    for(int i = 0; i < dim; i++) {
        printf(" x[ %d ] = %23.15e\n", i, x[i]);
    }*/
}