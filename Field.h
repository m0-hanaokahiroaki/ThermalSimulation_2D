#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
using namespace std;

struct Field {
protected:
	int Nx;//x�����̊i�q��
	int Ny;//y�����̊i�q��
	double dt;//����[�X�e�b�v]s]
	vector<vector<double>> alpha;//�̈�̔M�g�U��[m^2/s]
	const double alpha_w = 1.651 * 1e-7;//360K�̐��̔M�g�U��[m^2/s]
	const double T_h = 273.0 + 87.0;//���̉��x[K]
	const double epsilon = 1e-12;//������������̒l
	double h;//�i�q�� [m]
	vector<vector<double>> T;//���݂̉��x�� [K]
	vector<vector<double>> Tnew;//�X�V�p�̉��x�� [K]
	vector<vector<vector<double>>> T_list;//���x��̎��n���vector
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

//Field�̏������G�p�����[�^�̐ݒ�A���x��T,�X�V�p�̉��x��Tnew�̏�����Ԃ�ݒ肷��
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

//���xT_c, �M�g�U��alpha_c�̒����`���̂�ݒu����
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

//���E�����̐ݒ�
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

//��������
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

//Gauss-Seidel�@�Ŏ��ԃX�e�b�vn�̉��x������߂�
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

//���݂̉��x��T���擾����
vector<vector<double>> Field::get_T() {
	return T;
}

//���݂܂ł̉��x���vector���擾����
vector<vector<vector<double>>> Field::get_T_list() {
	return T_list;
}

//����(���ԃX�e�b�vn)�̉��x����t�@�C���ɏo�͂���
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