
#define _CRT_SECURE_NO_WARNINGS

#include <string>
#include <complex>
#include <cmath>

#define TATE 256
#define YOKO 256
#define SIZE 256*256
#define PI 3.141592653589


void error(char filename[]) {
	printf("File not exist:%s¥n", filename);
}


std::complex<double> * recursive_fft(std::complex<double> *, int);
std::complex<double> * recursive_ifft(std::complex<double> *, int);

int main(int argc, char* argv[]){
	FILE *fpin, *fpout;
	char filenamein[40], filenameout[40], c;
	int imagein[TATE][YOKO], imageout[TATE][YOKO];
	static int imageout1[SIZE];
	static double imagein1[SIZE];
	int i, j, k = 0;
	double a, b;

	sprintf(filenamein, "BARBARA.Y");
	if ((fpin = fopen(filenamein, "rb")) == NULL) { error(filenamein); }
	for (i = 0; i<TATE; i++) {
		for (j = 0; j<YOKO; j++) {
			imagein[i][j] = getc(fpin);
		}
	}
	fclose(fpin);

  /////////  水平方向のFFT   ////////////
  std::complex<double> horizontalFFT[TATE][YOKO];
	for (i = 0; i < TATE; i++) {
	  std::complex<double> arr[YOKO]; //初期化
		for (j = 0; j < YOKO; j++) {
			arr[j] = std::complex<double>(imagein[i][j], 0);
		}

		std::complex<double> *pArr;
		int len;

		len = sizeof(arr) / sizeof(arr[0]);
		pArr = arr;

		std::complex<double> *pfft = recursive_fft(pArr, len);
		for (int j = 0; j < YOKO; j++) {
			horizontalFFT[i][j] = *(pfft + j);                 /////  FFT  /////
		}
	}

  // 1回目の行列の転置 http://www.clg.niigata-u.ac.jp/~medimg/practice_medical_imaging/imgproc_scion/5fourier/index.htm
  std::complex<double> transposeFFT[YOKO][TATE];
	for (i = 0; i < YOKO; i++) {
		for (j = 0; j < TATE; j++) {
        transposeFFT[i][j] = horizontalFFT[j][i];
    }
  }

  /////////  縦方向のFFT   ////////////
  std::complex<double> verticalFFT[YOKO][TATE];
  for (i = 0; i < YOKO; i++) {
    std::complex<double> arr[TATE]; //初期化
    for (j = 0; j < TATE; j++) {
      arr[j] = transposeFFT[i][j];
    }

    std::complex<double> *pArr;
    int len;

    len = sizeof(arr) / sizeof(arr[0]);
    pArr = arr;

    std::complex<double> *pfft = recursive_fft(pArr, len);
    for (int j = 0; j < TATE; j++) {
      verticalFFT[i][j] = *(pfft + j);                 /////  FFT  /////
    }
  }

  // 2回目の行列の転置 http://www.clg.niigata-u.ac.jp/~medimg/practice_medical_imaging/imgproc_scion/5fourier/index.htm
  std::complex<double> outputFFT[TATE][YOKO];
  for (i = 0; i < TATE; i++) {
    for (j = 0; j < YOKO; j++) {
        outputFFT[i][j] = verticalFFT[j][i];
    }
  }

  // 0〜255に正規化するために最大値と最小値をまずは取得する
  double max_power = 0.0;
  double min_power = 0.0;
  for (i = 0; i<TATE; i++) {
     for (j = 0; j<YOKO; j++) {
          double real = outputFFT[i][j].real();
          double imag = outputFFT[i][j].imag();
          double power = log(real*real+imag*imag);
          if (power > max_power){
            max_power = power;
          }
          if (power < min_power){
            min_power = power;
          }
     }
  }

  sprintf(filenameout, "FFT_output.pgm");
  if ((fpout = fopen(filenameout, "wb")) == NULL) { error(filenameout); }
  fprintf(fpout, "P5¥n%d %d¥n255¥n", YOKO, TATE); // Windowsはこっち
  // fprintf(fpout, "P5\n%d %d\n255\n", YOKO, TATE); // Macはこっち
  for (i = 0; i<TATE; i++) {
     for (j = 0; j<YOKO; j++) {
          double real = outputFFT[i][j].real();
          double imag = outputFFT[i][j].imag();
          double power = log(real*real+imag*imag);
          power = (power-min_power)/(max_power-min_power)*255.0; //0〜255に正規化
          fprintf(fpout, "%c", (int)power);
     }
  }
  fclose(fpout);
 return 0;
}


std::complex<double> * recursive_fft(std::complex<double> *pArr, int len){
 std::complex<double>* arr = new std::complex<double>[len];
 // initilize array
 for (int i = 0; i < len; i++){
   arr[i] = *(pArr+i);
 }

 if (len == 1) {
   return arr;
   delete[] arr;
 }

 int half_len = len/2;
 std::complex<double>* even_arr = new std::complex<double>[half_len];
 std::complex<double>* odd_arr = new std::complex<double>[half_len];
 for (int i = 0; i < half_len; i++){
   even_arr[i] = arr[2*i];
   odd_arr[i] = arr[2*i+1];
 }

 // Exponent portion
 std::complex<double> expo = std::complex<double>(0, -2) * std::complex<double>(PI/len, 0);
 std::complex<double> omega_n = exp(expo);
 std::complex<double> omega = std::complex<double>(1, 0);

 std::complex<double> *py_even = recursive_fft(even_arr, half_len);   // *ついてるのはおまけでポインタではない
 std::complex<double> *py_odd = recursive_fft(odd_arr, half_len);
 std::complex<double>* y_even = new std::complex<double>[half_len];
 std::complex<double>* y_odd = new std::complex<double>[half_len];
 for (int i = 0; i < half_len; i++){
   y_even[i] = *(py_even+i);     // 配列 py[i]のデータをとってくる
   y_odd[i] = *(py_odd+i);
 }

 std::complex<double>* y = new std::complex<double>[len];
 for (int i = 0; i < half_len; i++){
   y[i] = y_even[i] + omega * y_odd[i];
   y[i+half_len] = y_even[i] - omega * y_odd[i];
   omega = omega * omega_n;
 }

 delete[] y_even;
 delete[] y_odd;
 delete[] arr;

 std::complex<double> *py = y;
 return py;
 delete[] y;
}
