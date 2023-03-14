/*
 * DOA.cpp
 *
 *  Created on: Mar 6, 2023
 *      Author: ashitd
 */


#include <iostream>
#include <armadillo>
#include <fstream>
#include <cstdlib>


#define n_signal 4.0 //number of signals
#define pi 3.141592653589793
#define N 200 // snapshots
#define M 10 // number of array elements (receiver )
#define lambda 150.0 // wavelength
#define snr 20 // Signal to noise ratio
using namespace  std;
using namespace arma;
wall_clock timer;

class Signal {
public:
	//mat doa = randn<mat>(1,n_signal,distr_param(-1,1));
	//doa = round(doa * 100.0);
	dmat doa = {{-75 , 0 , 40 , 60 }};
	dmat w = {{pi/4 , pi/3 , pi/2 , pi/5 }};
	int P = w.n_cols;
	cx_dmat D = zeros<cx_dmat>(P,M); // create zero matrix with no of signals x no of elemnts
	mat MSpace = regspace(0,M-1).t(); // Receiver array
	cx_dmat MSpace_complex = zeros<cx_mat>(MSpace.n_rows,MSpace.n_cols);
	dmat snapshot_matrix  = regspace(1,N).t();
};

int PlotOutput (mat angles , mat PMusic ) {
	ofstream DOAFile("DOA_data.csv");
	DOAFile << "Angle(theta/degree),spectrum function P(theta)/dB\n";
	for (int i = 0 ; i < (int) angles.n_cols ; ++i) {
		DOAFile << angles(0,i) << "," << PMusic(0,i) << "\n";
	}
	DOAFile.close();
	//Plot using python api
	system("python ./plot_doa.py ./DOA_data.csv");
	//
	return 0;
}

int main() {
	timer.tic();
	arma_rng::set_seed_random();
	cout << "test ok"<< endl;
	complex<double> ii(0,1); // get j
	double d = lambda / 2;
	Signal X; // get basic signal defination
	X.doa.print("DOA original:");
	X.doa = X.doa / 180 * pi; // redian conversion
	X.w = X.w.t(); // tranpose Frequency vector
	X.MSpace_complex.set_real(X.MSpace); // receiver array in complex form
	for (int k = 0 ; k < X.P ; ++k) {
		 X.D.row(k)	=  exp((- ii * 2.0 * pi * d * sin(X.doa(0,k)) / lambda ) *  X.MSpace_complex); // Assignment Matrix
	}
	X.D = X.D.t(); // transpose assignment matrix
	mat L = X.w * X.snapshot_matrix ; // distribute Frequency over snapshot space
	cx_dmat xx = 2.0 * exp(ii * L); // simulate signal
	cx_mat x = X.D * xx; // distribute signals over Receiver array
	// need to add Gaussian noise distribution
	cx_mat R = x * x.t() ; // Data covarience matrix
	vec V;
	cx_mat NM;
    eig_sym(V,NM,R); // Find the eigenvalues and eigenvectors of R
    cx_mat NN = NM.cols(0,M-X.P-1); // Estimate noise subspace
    dmat theta = regspace(-90,0.5,90).t();
    cx_mat SS = zeros<cx_mat>(1,M);
    mat PMusic = zeros<mat>(1,(int) theta.n_cols);
    for (int i = 0 ; i < (int) theta.n_cols ; ++i ) {
    	SS = zeros<cx_mat>(1,M);
    	for (int j = 0 ; j < M ; ++j ) {
    		//cout << i << "---" << j << endl;
    		SS(0,j) = exp(-ii * 2.0 * (double) j * pi * d * sin(theta(0,i)/180*pi) / lambda);
    	}
    cx_mat PP = SS * NN * NN.t() * SS.t();
    PP = 1.0 / PP;
    mat PP_abs = abs(PP);
    PMusic(0,i) = PP_abs(0,0);
    }
    PMusic = 10 * log10(PMusic / PMusic.max());
    double time_taken = timer.toc();
    cout << "time taken " << time_taken << " seconds "<< endl;
    PlotOutput(theta,PMusic);
}
