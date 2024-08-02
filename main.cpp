#include <iostream>
#include <vector>
#include <iomanip>
#include <chrono>
#include "Field.h"
using namespace std;

const int Nx = 101;//x方向の格子数
const int Ny = 101;//y方向の格子数
const double dt = 0.1;//時間ステップ幅[s]
const double h = 0.003;//格子幅[m]
const double alpha_c = 2.27 * 1e-5;//物体の熱拡散率[m^2/s]
const double T_c = 273.0 + 15.0;//物体の温度[K]
const double T_end = 273.0 + 50.0;// 中心温度の目標温度[K]

int n_max = 10000;//最大時間ステップ数
const bool output = true;//結果をファイルに出力するかどうか

int main() {

	//計算時間測定開始
	auto t0 = chrono::high_resolution_clock::now();

	//Fieldの宣言
	Field field;
	
	//Fieldの初期化
	field.init(Nx, Ny, dt, h, alpha_c, T_c);

	//出力する小数点を一桁に固定
	cout << fixed << setprecision(1);
	//0ステップ目(初期状態)の温度場をファイルに出力する
	if (output)field.write_file(0);

	//それぞれの時間ステップごとにGauss-Seidel法で温度場を求める
	for (int n = 0; n < n_max; n++) {
		
		//nステップ目の温度場を求める
		field.GaussSeidel(n);

		//nステップ目の温度場をファイルに出力する
		if(output)field.write_file(n+1);
		
		//領域の中心温度がT_end以上かどうか(T_end以上ならば時間ステップ終了、T_end未満ならば時間ステップを進める)
		if (field.get_T()[Nx / 2][Ny / 2] >= T_end) n_max = n + 1;

		//100時間ステップごとに中心温度を確認
		if (n % 100 == 0) cout << "t = " << n*dt << "s : 中心温度 = " << field.get_T()[Nx / 2][Ny / 2] - 273.0 << "℃" << endl;
	}
	cout << endl;
	if (output)cout << "output_0.txt～outputfile_" << n_max << ".txt were successfully saved." << endl;

	//中心温度と時間tの確認
	cout << "中心温度 = " << field.get_T()[Nx / 2][Ny / 2] - 273.0 << "℃ は t = " << n_max * dt << "s で達成されました。" << endl;

	//計算時間測定終了
	auto t1 = chrono::high_resolution_clock::now();
	cout << "計算時間は " << double(chrono::duration_cast<chrono::milliseconds>(t1 - t0).count()) / 1000.0 << " s";

	return 0;
}