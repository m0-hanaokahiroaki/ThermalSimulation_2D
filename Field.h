#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
using namespace std;

struct Field {
protected:
	int Nx;//x方向の格子数
	int Ny;//y方向の格子数
	double dt;//時間[ステップ]s]
	vector<vector<double>> alpha;//領域の熱拡散率[m^2/s]
	const double alpha_w = 1.651 * 1e-7;//360Kの水の熱拡散率[m^2/s]
	const double T_h = 273.0 + 87.0;//水の温度[K]
	const double epsilon = 1e-12;//収束判定条件の値
	double h;//格子幅 [m]
	vector<vector<double>> T;//現在の温度場 [K]
	vector<vector<double>> Tnew;//更新用の温度場 [K]
	vector<vector<vector<double>>> T_list;//温度場の時系列のvector
public:
	void init(int newNx, int newNy, double newdt, double newh, double alpha_c, double T_c);
	void setObject(double alpha_c, double T_c);
	void BoundaryCondition();
	bool IsConverged();
	void GaussSeidel(int n);
	vector<vector<double>> get_T();
	vector<vector<vector<double>>> get_T_list();
	void write_file(int n);
};

//Fieldの初期化；パラメータの設定、温度場T,更新用の温度場Tnewの初期状態を設定する
void Field::init(int newNx, int newNy, double newdt, double newh, double alpha_c, double T_c) {
	Nx = newNx;
	Ny = newNy;
	dt = newdt;
	h = newh;
	T.resize(Nx, vector<double>(Ny));
	alpha.resize(Nx, vector<double>(Ny));
	setObject(alpha_c, T_c);
	BoundaryCondition();
	Tnew = T;
	T_list.push_back(T);
}

//温度T_c, 熱拡散率alpha_cの長方形物体を設置する
void Field::setObject(double alpha_c, double T_c) {
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			if (i == 0 || i == Nx - 1 || j == 0 || j == Ny - 1) {
				alpha[i][j] = alpha_w;
				T[i][j] = T_h;
			}
			else {
				alpha[i][j] = alpha_c;
				T[i][j] = T_c;
			}
		}
	}
}

//境界条件の設定
void Field::BoundaryCondition() {
	for (int i = 0; i < Nx; i++) {
		T[i][0] = T_h;
		T[i][Ny - 1] = T_h;
	}
	for (int j = 0; j < Ny; j++) {
		T[0][j] = T_h;
		T[Nx - 1][j] = T_h;
	}
}

//収束判定
bool Field::IsConverged() {
	double r_sum = 0.0;
	double Tnew_sum = 0.0;
	for (int i = 1; i < Nx - 1; i++) {
		for (int j = 1; j < Ny - 1; j++) {
			r_sum += abs(Tnew[i][j] - T[i][j]);
			Tnew_sum += abs(Tnew[i][j]);
		}
	}
	if (r_sum / Tnew_sum < epsilon) {
		return true;
	}
	else {
		return false;
	}
}

//Gauss-Seidel法で時間ステップnの温度場を求める
void Field::GaussSeidel(int n) {
	while (true) {
		BoundaryCondition();
		T = Tnew;
		for (int i = 1; i < Nx - 1; i++) {
			for (int j = 1; j < Ny - 1; j++) {
				Tnew[i][j] = (alpha[i][j] * dt / (2 * h * h)) * (T[i + 1][j] + Tnew[i - 1][j] + T[i][j + 1] + Tnew[i][j - 1] - 4 * T[i][j]
					+ T_list[n][i + 1][j] + T_list[n][i - 1][j] + T_list[n][i][j + 1] + T_list[n][i][j - 1] - 4 * T_list[n][i][j]) + T_list[n][i][j];
			}
		}
		if (IsConverged())break;
	}
	T = Tnew;
	T_list.push_back(T);
}

//現在の温度場Tを取得する
vector<vector<double>> Field::get_T() {
	return T;
}

//現在までの温度場のvectorを取得する
vector<vector<vector<double>>> Field::get_T_list() {
	return T_list;
}

//現在(時間ステップn)の温度場をファイルに出力する
void Field::write_file(int n) {

	string fileName = "output_"+to_string(n)+".txt";
	ofstream outfile(fileName);

	outfile << fixed << setprecision(2);
	for (auto i = T.begin(); i != T.end(); ++i) {
		for (auto ij = i->begin(); ij != i->end(); ++ij) {
			outfile << *ij - 273.15 << " ";
		}
		outfile << endl;
	}
	outfile.close();

}