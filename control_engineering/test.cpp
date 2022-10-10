// control engineering

#include <cmath>
#include <GL/freeglut.h>
#include <fstream>
#include <vector>
#include <stdio.h>

namespace{
	using myvec = std::vector<double>;

	static auto constexpr VARIABLE = 4;



	static auto constexpr TIME_PERIOD = 1; // タイマ周期 ミリ秒
	static auto constexpr DT = 0.001;  // 数値積分の時間刻み 秒

	static auto constexpr c = 0.1; // 振子の粘性抵抗
	static auto constexpr F = 0.1; // 台車の粘性抵抗
	static auto constexpr l = 1.0; // 振子の長さ
	static auto constexpr m = 1.0; // 振子の質量
	static auto constexpr M = 1.0; // 台車の質量
	static auto constexpr J = 1.0; // 振子の慣性モーメント
	static auto constexpr g = 9.8; // 重力加速度

	// 変数
	// auto x0 = 0.0;
	// auto xdot0 = 0.1;
	// auto theta0 = 0.0;
	// auto thetadot0 = 0.0;

	// myvec q{x0, xdot0, theta0, thetadot0};

	double x[4]  = {0.0, 0.1, 0.0, 0.0}; // 状態変数
	double dx[4] = {0.0, 0.0, 0.0, 0.0}; // 状態変数の微分
	double u     = 0.0; // 制御入力 台車の駆動力
	double t     = 0.0; // 現在時刻

	// ウィンドウサイズ
	int window_wiDTh;
	int window_height;

	// ファイル
	FILE* file = NULL;

	// タイマコールバック関数
	void timer_func(int value);

	// キーボードコールバック関数
	void keyboard_func(unsigned char ch, int x, int y);

	// ウィンドウサイズ変更コールバック関数
	void reshape_func(int wiDTh, int height);

	// 画面描画コールバック関数
	void display_func();

	// ウィンドウが閉じられたときに呼ばれるコールバック関数
	void close_func();
}

int main(int argc, char** argv){
	// OpenGLを初期化
	glutInit(&argc, argv);

	// ウィンドウを作成
	glutCreateWindow("pendulum");

	// コールバック関数を登録
	glutTimerFunc(TIME_PERIOD, timer_func, 0);
	glutKeyboardFunc(keyboard_func);
	glutReshapeFunc(reshape_func);
	glutDisplayFunc(display_func);

	// ファイルを開く
	file = fopen("main_output.txt", "w");
	
	// メインループ（無限ループ）
	glutMainLoop();

	return 0;
}

namespace{
	void timer_func(int value){
		// 状態変数から各成分を取得
		double z      = x[0];
		double theta  = x[1];
		double dz     = x[2];
		double DTheta = x[3];

		// ここで制御入力を決定
		double u = 0.0;

		// ここに運動方程式を記述する
		// double A[2][2];
		// double b[2];

		double A[2][2] = { {M + m, m * l * cos(theta)}, {m * l * cos(theta), J + m * pow(l,2)} };
 		double B[2] = {m * pow(l,2) * pow(DTheta,2) * sin(theta) - F * dz, -c * DTheta + m * g * l * sin(theta)};

		// A[0][0] = M + m;
		// A[0][1] = m * l * cos(theta);
		// A[1][0] = m * l * cos(theta);
		// A[1][1] = J + m * pow(l,2);

		double tau[2] = {u,0};

		double a_inv = 1 / (A[0][0] * A[1][1] - A[0][1] * A[1][0]);

		double A_inv[2][2] = { {A[1][1] * a_inv, -A[0][1] * a_inv }, {-A[1][0] * a_inv, A[0][0] * a_inv} };

		// 台車加速度
		double ddz = A_inv[0][0] * (B[0] + tau[0]) + A_inv[0][1] * (B[1] + tau[1]);
		// 振り子角速度
		double dDTheta = A_inv[1][0] * (B[0] + tau[0]) + A_inv[1][1] * (B[1] + tau[1]);

		// 状態変数の時間微分
		dx[0] = dz;
		dx[1] = DTheta;
		dx[2] = ddz;
		dx[3] = dDTheta;
	
		// 数値積分（オイラー則）
		for(int i = 0; i < 4; i++){
			x[i] += dx[i] * DT;
		}

		// 時刻の更新
		t += DT;

		// ファイルに出力
		fprintf(file, "%f, %f, %f, %f, %f\n", t, x[0], x[1], x[2], x[3]);

		// 再描画要求
		glutPostRedisplay();

		// タイマコールバックを再登録
		glutTimerFunc(TIME_PERIOD, timer_func, 0);
	}

	void keyboard_func(unsigned char ch, int x, int y){

	}

	void reshape_func(int wiDTh, int height){
		window_wiDTh  = wiDTh;
		window_height = height;

		// 描画範囲を設定
		glViewport(0, 0, wiDTh, height);
	}

	void display_func(){
		double w = (double)window_wiDTh ;
		double h = (double)window_height;

		// 描画処理のお膳立て
		glClear(GL_COLOR_BUFFER_BIT);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(-w/2.0, w/2.0, -h/2.0, h/2.0, -1.0, 1.0);

		// 座標軸を描画
		glBegin(GL_LINES);
		glColor3f(1.0f, 1.0f, 1.0f);
		glVertex2d(-w/2.0, 0.0);
		glVertex2d( w/2.0, 0.0);
		glVertex2d(0.0, -h/2.0);
		glVertex2d(0.0,  h/2.0);
		glEnd();

		auto const scale = 100.0;

		// 振り子を描画
		glBegin(GL_LINES);
		glColor3f(1.0f, 1.0f, 1.0f);
		glVertex2d(scale * x[0], 0.0);
		glVertex2d(scale * (x[0] + l * sin(x[1])), scale * l * cos(x[1]));
		glEnd();

		glutSwapBuffers();
	}

	void close_func(){
		// ファイルを閉じる
		fclose(file);
	}
}
