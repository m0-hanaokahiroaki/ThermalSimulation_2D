#include <iostream>
#include <vector>
#include <iomanip>
#include <chrono>
#include "Field.h"
using namespace std;

const int Nx = 101;//x�����̊i�q��
const int Ny = 101;//y�����̊i�q��
const double dt = 0.1;//���ԃX�e�b�v��[s]
const double h = 0.003;//�i�q��[m]
const double alpha_c = 2.27 * 1e-5;//���̂̔M�g�U��[m^2/s]
const double T_c = 273.0 + 15.0;//���̂̉��x[K]
const double T_end = 273.0 + 50.0;// ���S���x�̖ڕW���x[K]

int n_max = 10000;//�ő厞�ԃX�e�b�v��
const bool output = true;//���ʂ��t�@�C���ɏo�͂��邩�ǂ���

int main() {

	//�v�Z���ԑ���J�n
	auto t0 = chrono::high_resolution_clock::now();

	//Field�̐錾
	Field field;
	
	//Field�̏�����
	field.init(Nx, Ny, dt, h, alpha_c, T_c);

	//�o�͂��鏬���_���ꌅ�ɌŒ�
	cout << fixed << setprecision(1);
	//0�X�e�b�v��(�������)�̉��x����t�@�C���ɏo�͂���
	if (output)field.write_file(0);

	//���ꂼ��̎��ԃX�e�b�v���Ƃ�Gauss-Seidel�@�ŉ��x������߂�
	for (int n = 0; n < n_max; n++) {
		
		//n�X�e�b�v�ڂ̉��x������߂�
		field.GaussSeidel(n);

		//n�X�e�b�v�ڂ̉��x����t�@�C���ɏo�͂���
		if(output)field.write_file(n+1);
		
		//�̈�̒��S���x��T_end�ȏォ�ǂ���(T_end�ȏ�Ȃ�Ύ��ԃX�e�b�v�I���AT_end�����Ȃ�Ύ��ԃX�e�b�v��i�߂�)
		if (field.get_T()[Nx / 2][Ny / 2] >= T_end) n_max = n + 1;

		//100���ԃX�e�b�v���Ƃɒ��S���x���m�F
		if (n % 100 == 0) cout << "t = " << n*dt << "s : ���S���x = " << field.get_T()[Nx / 2][Ny / 2] - 273.0 << "��" << endl;
	}
	cout << endl;
	if (output)cout << "output_0.txt�`outputfile_" << n_max << ".txt were successfully saved." << endl;

	//���S���x�Ǝ���t�̊m�F
	cout << "���S���x = " << field.get_T()[Nx / 2][Ny / 2] - 273.0 << "�� �� t = " << n_max * dt << "s �ŒB������܂����B" << endl;

	//�v�Z���ԑ���I��
	auto t1 = chrono::high_resolution_clock::now();
	cout << "�v�Z���Ԃ� " << double(chrono::duration_cast<chrono::milliseconds>(t1 - t0).count()) / 1000.0 << " s";

	return 0;
}