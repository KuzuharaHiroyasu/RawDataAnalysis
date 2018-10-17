// RawDataAnalysis.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//

#include	"stdafx.h"
#include	"hanning_window.h"
#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	<time.h>
#include	<direct.h>
#include	<string.h>

/*==============================================================================*/
/*	define定義																	*/
/*==============================================================================*/
#define _SAMPLING_FREQUENCY			0.15625f
#define _PULSE_PARAM_SNPK_COEF		0.05f
#define _PULSE_PARAM_START			6
#define _PULSE_PARAM_END			26
#define BUF_SIZE					128
#define	pi							((double)3.1415926535)

/*==============================================================================*/
/*	グローバル変数																*/
/*==============================================================================*/
int	    temp_int_buf0[BUF_SIZE];
double	temp_dbl_buf0[BUF_SIZE];
double	temp_dbl_buf1[BUF_SIZE];
double	temp_dbl_buf2[BUF_SIZE];
double	pdata[BUF_SIZE];
double  ifft_[BUF_SIZE];
double  fft_[BUF_SIZE];
char    path_[17] = "./";

int len = 128;

double ac_avg_clr;
double ac_avg_inf;
double dc_avg_clr;
double dc_avg_inf;

/*==============================================================================*/
/*	プロトタイプ宣言															*/
/*==============================================================================*/
void fft(const double in[], int N, double ar[], double ai[], double p[]);
void ifft(double	ar[], double	ai[], int N, double	ot[]);
void peak_vallay(double	in[], int ot[], int size, int width, double	th, int	peak);
int	 peak_modify(double in_data[], int in_res[], double ot_data[], double ot_hz[], int size, double delta);	/* ☆ */
void MOD(double PM_M, double PM_0, double PM_P, int m, double DELTA, double *FI, double *PP);
void acdc_average(const double* pdata, double* ar1, double* ai1, double* p3, double pulse, int len, int no);
void debug_out(char *f, const double d[], int size, const char* ppath, int no);
void debug_out_int(char *f, const int d[], int size, const char* ppath, int no);

/*==============================================================================*/
/*	main																		*/
/*==============================================================================*/
int main()
{
	FILE *fp;
	const double coeff = _PULSE_PARAM_SNPK_COEF;
	const int start = _PULSE_PARAM_START;
	const int end = _PULSE_PARAM_END;
	struct tm timeptr;
	char folder[15] = { '\0' };

	double* ph11 = NULL;
	double* p2 = NULL;
	double* p3 = NULL;
	double*	p4 = NULL;
	double* ar2 = NULL;
	double* ai2 = NULL;
	double*	f1 = NULL;
	int* peak1 = NULL;

	time_t timer = 0;
	double max = 0;
	int i = 0;
	int ii = 0;
	int no = 0;

	if (!fopen_s(&fp, "./raw(1).txt", "r")) {
		no = 1;
	}
	else if (!fopen_s(&fp, "./raw(2).txt", "r")) {
		no = 2;
	}
	else {
		printf("ファイルを開くことが出来ませんでした。\n");
		return 1;
	}

	timer = time(NULL);
	if (localtime_s(&timeptr, &timer))
	{
		fclose(fp);
		return 1;
	}

	strftime(folder, 15, "%Y%m%d%H%M%S", &timeptr);

	if (!strcat_s(path_, sizeof path_, folder)) {
		_mkdir(path_);

		for (i = 0; i < len; i++) {
			fscanf_s(fp, "%lf", &(pdata[i]));     /*  1行読む  */
		}

		/*- DC成分除去 -----------------------------------------------------------------*/
		ph11 = &temp_dbl_buf1[0];
		double min = pdata[0];
		for (int ii = 1; ii < len; ++ii) {
			if (min > pdata[ii]) {
				min = pdata[ii];
			}
		}
		for (int ii = 0; ii < len; ++ii) {
			ph11[ii] = (pdata[ii] - min);			//DC抜きデータ
		}
		//		debug_out("raw", pdata, len, path_, no);	//生データ出力
		debug_out("dc", ph11, len, path_, no);		//DC抜きデータ出力

													/*- fft --------------------------------------------------------------------*/
		double	ar1[BUF_SIZE];		//実数部
		double	ai1[BUF_SIZE];		//虚数部

		for (int ii = 0; ii < len; ++ii) {
			ph11[ii] = pdata[ii] * hanning_window[ii];		//演算データ
		}
		fft(ph11, len, ar1, ai1, fft_);
		debug_out("fft", fft_, len, path_, no);		//FFT演算結果出力
		debug_out("ar1", ar1, len, path_, no);		//実数部出力
		debug_out("ai1", ai1, len, path_, no);		//虚数部出力

													/*- 不要なデータをマスクする------------------------------------------------*/
													/*- パワーを求める----------------------------------------------------------*/
		p2 = &temp_dbl_buf1[0];
		// startまではall 0
		for (ii = 0; ii < start; ++ii) {
			p2[ii] = 0;
		}
		// ２乗根処理(振幅スペクトル)
		for (ii = start; ii < end; ++ii) {
			p2[ii] = sqrt((sqrt(ar1[ii] * ar1[ii] + ai1[ii] * ai1[ii])) / len);
		}
		for (ii = end; ii < len; ++ii) {
			p2[ii] = 0;
		}
		debug_out("p2", p2, len, path_, no);		//帯域抽出出力

													/*- ピーク検出 --------------------------------------------------------------*/
		peak1 = &temp_int_buf0[0];	//ピーク判定結果
		p4 = &temp_dbl_buf0[0];		//パワースペクトル
		f1 = &temp_dbl_buf2[0];		//周波数

									//ピーク判定
		peak_vallay(p2, peak1, len, 3, 0.1, 1);
		//パワースペクトル、周波数出力
		peak_modify(p2, peak1, p4, f1, len, 0.1);

		/*- 逆FFT --------------------------------------------------------------*/
		ar2 = &temp_dbl_buf0[0];	//実数部
		ai2 = &temp_dbl_buf2[0];	//虚数部

		memcpy(&ar2[0], &p2[0], sizeof(double)*len);	//実数部に帯域抽出値コピー
		memset(&ai2[0], 0x00, len);		//虚数部に０セット

		p3 = &ifft_[0];
		ifft(ar2, ai2, len, p3);
		debug_out("ifft", p3, len, path_, no);			//逆FFT演算結果出力
		debug_out("ifft_ar2", ar2, len, path_, no);		//実数部出力
		debug_out("ifft_ai2", ai2, len, path_, no);		//虚数部出力

														/*- 最大値との比を計算 --------------------------------------------------------------*/
		max = 0;
		for (ii = start; ii < end; ++ii) {
			if (max < p3[ii]) {
				max = p3[ii];
			}
		}
		for (ii = 0; ii < len; ++ii) {
			p3[ii] /= max;
		}
		debug_out("p3hi", p3, len, path_, no);		//最大値との比出力

													/*- ピーク検出 --------------------------------------------------------------*/
		peak1 = &temp_int_buf0[0];
		p4 = &temp_dbl_buf1[0];
		f1 = &temp_dbl_buf2[0];
		peak_vallay(p3, peak1, len, 3, 0.1, 1);
		peak_modify(p3, peak1, p4, f1, len, 1);
		debug_out("peak", f1, len, path_, no);

		/*- HR --------------------------------------------------------------*/
		double pulse = (int)(60 / (f1[0] * coeff));
		debug_out("snpk", &pulse, 1, path_, no);

		acdc_average(pdata, ar1, ai1, p3, pulse, len, no);	//AC、DCの平均値算出
	}

	fclose(fp);
}

/*==============================================================================*/
/*	fft																			*/
/*==============================================================================*/
void fft
(
	const double	in[],	/* IN：実数部										*/
	int		N,				/* IN：要素数										*/
	double	ar[],			/* OT：実数部										*/
	double	ai[],			/* OT：虚数部										*/
	double	p[]				/* OT：パワー										*/
)
{
	double	ReF = 0.0;
	double	ImF = 0.0;

	// 実数部分と虚数部分に分けてフーリエ変換
	for (int n = 0; n<N; n++)
	{
		ReF = ImF = 0.0;
		for (int k = 0; k<N; k++)
		{
			ReF += in[k] * cos(2 * pi*k*n / N);
			ImF += -in[k] * sin(2 * pi*k*n / N);
		}
		ar[n] = ReF;
		ai[n] = ImF;
		p[n] = 0;
	}
}

/*==============================================================================*/
/*	ifft																		*/
/*==============================================================================*/
void ifft
(
	double	ar[],	/* IN：実数部												*/
	double	ai[],	/* IN：虚数部												*/
	int		N,		/* IN：要素数												*/
	double	ot[] 	/* OT：逆生成データ											*/
)
{
	int		n;
	int		k;
	double	Ref;

	for (k = 0; k<N; k++)
	{
		Ref = 0.0;
		for (n = 0; n<N; n++)
		{
			Ref += (ar[n] * cos(2 * pi*k*n / N) - ai[n] * sin(2 * pi*k*n / N));
		}
		Ref /= ((double)N);
		ot[k] = Ref;
	}
}

/*==============================================================================*/
/*	peak_vallay																	*/
/*==============================================================================*/
void peak_vallay
(
	double	in[],		/* IN：入力信号の格納されたバッファ						*/
	int		ot[],		/* OT：判定結果（0=NO、1=YES)							*/
	int		size,		/* IN：入力信号・出力信号のバッファのサイズ				*/
	int		width,		/* IN：検査幅、奇数、最小は３指定						*/
	double	th,			/* IN：検査閾値											*/
	int		peak		/* IN：1=極大値抽出、else=極小処理抽出					*/
)
{
	int		w = (width / 2);

	/*--------------------------------------------------------------------------*/
	for (int i = 0; i < size; i++)
	{
		ot[i] = 0;
	}

	// width = 3 -> w = 1 前提のコード
	int loop = size - 1;
	for (int ii = 1; ii<loop; ++ii) {
		if (peak == 1) {
			if ((in[ii - 1] < in[ii]) && (in[ii] > in[ii + 1])) {
				if (in[ii] > th) {
					ot[ii] = 1;
				}
			}
		}
		else {
			if ((in[ii - 1] > in[ii]) && (in[ii] < in[ii + 1])) {
				if (in[ii] > th) {
					ot[ii] = 1;
				}
			}
		}
	}
}

/*==============================================================================*/
/*	peak_modify	極大値補正														*/
/*==============================================================================*/
int peak_modify																	/* ☆ */
(
	double	in_data[],	/* IN：入力信号の格納されたバッファ						*/
	int		in_res[],	/* IN：判定結果（0=NO、1=YES)							*/
	double	ot_data[],	/* OT：出力信号パワースペクトル							*/
	double	ot_hz[],	/* OT：出力信号周波数									*/
	int		size,		/* IN：入力信号・出力信号のバッファのサイズ				*/
	double	delta		/* IN：FFTサンプリング周波数							*/
)
{
	int pos = 0;
	/*--------------------------------------------------------------------------*/
	for (int i = 0; i < size; i++)
	{
		ot_data[i] = 0.0;
		ot_hz[i] = 0.0;
	}

	/*--------------------------------------------------------------------------*/
	for (int i = 1; i < (size - 1); i++)
	{
		if (in_res[i] == 1)
		{
			double PM_M = in_data[i - 1];
			double PM_0 = in_data[i];
			double PM_P = in_data[i + 1];
			double FI = 0.0;
			double PP = 0.0;

			MOD(PM_M, PM_0, PM_P, i, delta, &FI, &PP);

			ot_data[pos] = PP;
			ot_hz[pos] = FI;
			pos += 1;
		}
	}

	return pos;
}

/*==============================================================================*/
/*	modify		極大値補正														*/
/*==============================================================================*/
/*==============================================================================*
[入力]PM_M  : ローカルピーク周波数のひとつ前の周波数（パワースペクトル）
[入力]PM_0  : ローカルピーク周波数（パワースペクトル）
[入力]PM_P  : ローカルピーク周波数のひとつ後の周波数（パワースペクトル）
[入力]DELTA : サンプリング周波数（１単位）＝定数
[入力]m     : ローカルピーク周波数（パワースペクトル）の位置＝データの配列番号に等しい
[出力]FI    : 算出ピーク値の周波数
[出力]PP    : 算出ピーク値のパワースペクトル
*==============================================================================*/
static	void	MOD																	/* ☆ */
(
	double PM_M,
	double PM_0,
	double PM_P,
	int    m,
	double DELTA,
	double *FI,
	double *PP
)
{
	double a = 0.0;
	double b = 0.0;
	double c = 0.0;

	*FI = 0.0;
	*PP = 0.0;

	a = ((PM_P + PM_M) / 2.0) - PM_0;			/* (3.41)参照				*/
	b = ((PM_P - PM_M) / 2.0);				/* (3.41)参照				*/
	c = PM_0;										/* (3.41)参照				*/

	*FI = ((-b) / (2.0 * a) + m) * DELTA;		/* (3.44)参照				*/
	*PP = c - ((b * b) / (4.0 * a));			/* (3.45)参照				*/
}

/*==============================================================================*/
/*	acdc_average																*/
/*==============================================================================*/
static void acdc_average(const double* pdata, double* ar1, double* ai1, double* p3, double pulse, int len, int no)
{
	double ac_avg = 0;
	double dc_avg = 0;
	double pos_center = ((pulse / 60) / _SAMPLING_FREQUENCY);	//センター位置
	int pos_center_down;
	int ii;

	// センターが0～127の範囲にない場合は0に丸める(不正アクセスエラー対策)
	if (pos_center < 0 || pos_center > len - 1) {
		pos_center = 0;
	}

	pos_center_down = (int)pos_center;

	// pos_center_downまでは 0
	for (int ii = 0; ii<pos_center_down; ++ii) {
		ar1[ii] = 0;
		ai1[ii] = 0;
	}

	if (pos_center_down == pos_center)
	{
		pos_center_down++;
	}
	else {
		pos_center_down = pos_center_down + 2;
	}
	// pos_center_down以降も 0
	for (ii = pos_center_down; ii<len; ++ii) {
		ar1[ii] = 0;
		ai1[ii] = 0;
	}
	debug_out("avg_ar1", ar1, len, path_, no);
	debug_out("avg_ai1", ai1, len, path_, no);
	// 逆FFT
	ifft(ar1, ai1, len, p3);
	debug_out("avg_ifft", p3, len, path_, no);

	// 演算結果を絶対値にする
	for (ii = 0; ii < len; ++ii) {
		p3[ii] = fabs(p3[ii]);
	}

	// 平均算出
	for (ii = 0; ii < len; ii++) {
		ac_avg += p3[ii];
		dc_avg += pdata[ii];
	}
	ac_avg /= len;
	dc_avg /= len;

	if (no == 1) {
		ac_avg_clr = ac_avg;	//赤色のAC平均値
		dc_avg_clr = dc_avg;	//赤色のDC平均値
	}
	else {
		ac_avg_inf = ac_avg;	//赤外のAC平均値
		dc_avg_inf = dc_avg;	//赤外のDC平均値
	}
	debug_out("ac_avg", &ac_avg, 1, path_, no);
	debug_out("dc_avg", &dc_avg, 1, path_, no);
}

/*==============================================================================*/
/*	debug_out																	*/
/*==============================================================================*/
void	debug_out(char *f, const double d[], int size, const char* ppath, int no)
{
	FILE		*fp;
	char		b[1024];

	sprintf_s(b, sizeof(b), "%s\\%s(%d).txt", ppath, f, no);


	if (fopen_s(&fp, b, "w"))
	{
		printf("file open error [debug_out]\n");
		exit(0);
	}

	for (int i = 0; i < size; i++)
	{
		fprintf(fp, "%lf\n", d[i]);
	}
	fclose(fp);
}

/*==============================================================================*/
/*	debug_out_int																*/
/*==============================================================================*/
void	debug_out_int(char *f, const int d[], int size, const char* ppath, int no)
{
	FILE		*fp;
	char		b[1024];

	sprintf_s(b, sizeof(b), "%s\\%s(%d).txt", ppath, f, no);


	if (fopen_s(&fp, b, "w"))
	{
		printf("file open error [debug_out]\n");
		exit(0);
	}

	for (int i = 0; i < size; i++)
	{
		fprintf(fp, "%d\n", d[i]);
	}
	fclose(fp);
}
/* EOF */
