
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

	std::complex<double> arr[YOKO];

	for (i = 0; i < TATE; i++) {            /////////  縦方向   ////////////
		for (j = 0; j < YOKO; j++) {
			arr[j] = std::complex<double>(imagein[i][j], 0);
		}

		std::complex<double> *pArr;
		int len;

		len = sizeof(arr) / sizeof(arr[0]);
		pArr = arr;

		std::complex<double> *pfft = recursive_fft(pArr, len);
		std::complex<double>* fft = new std::complex<double>[len];
		for (int j = 0; j < len; j++) {
			fft[j] = *(pfft + j);                 /////  FFT  /////
		}

		for (int j = 0; j < len; j++) {
			a = pow((double)fft[j].real(), 2.0);
			b = pow((double)fft[j].imag(), 2.0);
			imageout[i][j] = (int)log(sqrt(a + b));

			printf("imageout[%d][%d] = %d¥n", i,j,imageout[i][j]);

			if (imageout[i][j] < 0) {
				imageout[i][j] = 0;
			}
			if (imageout[i][j] > 255) {
				imageout[i][j] = 255;
			}
		}
	}

  sprintf(filenameout, "FFT_output.pgm");
  if ((fpout = fopen(filenameout, "wb")) == NULL) { error(filenameout); }
  if ((fpout = fopen(filenameout, "wb")) == NULL) {}
  fprintf(fpout, "P5¥n%d %d¥n255¥n", YOKO, TATE);
  for (i = 0; i<TATE; i++) {
     for (j = 0; j<YOKO; j++) {
         fprintf(fpout, "%c", imageout[i][j]);
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
