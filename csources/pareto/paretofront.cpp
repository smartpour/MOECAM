// paretofront.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <random>
#include <chrono>
#include <string>
#include <list>
#include <set>

	using namespace std;


	//	ofstream outputFile;

	ofstream outputFile1;
	ofstream outputFile2;
	ofstream outputFile3;


	int main(int argc, char **argv)

	{

		string inputString;
		string outputString;

		long int NData = 125;
		int dim = 2;
		int strict = 0;
		int cols = 0;

		if (argc != 7)
		{
			cout << "Usage: " << argv[0] << " inputfile outputfile dim ndata strict\n";
			system("pause");
			exit(0);
		}

		inputString = string(argv[1]);
		outputString = string(argv[2]);
		dim = atoi(argv[3]);
		NData = atol(argv[4]);
		strict = atol(argv[5]);
		cols = atol(argv[6]);

		ifstream inputFile(inputString);
		outputFile1.open(outputString);
		outputFile2.open(outputString+".label");

		outputFile3.open(outputString + ".test");

		cout << dim << " " << NData << endl;


		vector<float> temp(dim);

		vector<float> Data;

		vector<float> Pareto;
		set<long int> Paretoindex;
		list<long int> Erase;
		vector<long int> x;
		float r;
		long int t;

		for (long int j = 0; j < NData;j++) {  //GB1
			t = j;
			if (cols == 1) inputFile >> t;
			if (cols == 2) {inputFile >> t; inputFile >> t;}
			x.push_back(t);
		for (int i = 0; i < dim; i++) {
			inputFile >> r;
			Data.push_back(r);
			//cout << r << " ";
		}
	}
//		cout << Data.size()<<" "<<x.size()<<endl;
//		for (int j = 0; j < NData; j++) cout << x[j] << " ";
//		for (int k = 0;k < NData;k++) {
//			for (int j = 0; j < dim; j++)  outputFile1 << Data[dim*x[k] + j] << " ";
//			outputFile1 << endl;
//		}


		inputFile.close();

		Paretoindex.insert(0);

		long int d, nd;

		for (long int i = 1;i < NData;i++) {
			d = 0;
			Erase.clear();

			if (i % 10000 == 0) cout << i<< endl;
			for (auto it = Paretoindex.begin(); it != Paretoindex.end(); it++) {
				nd = 0;
				if (strict == 1) { // strong pareto
					for (int j = 0;j < dim;j++) {
						if (Data[i*dim + j] >= Data[*it  *dim + j]) //not dominated    GB3 was >-
						{
							nd++;
						}
					}
				}
				else { // weak pareto
					for (int j = 0;j < dim;j++) {
						if (Data[i*dim + j] > Data[*it  *dim + j]) //not dominated    GB3 was >-
						{
							nd++;
						}
					}
				}
				if(nd==0) d = 1; // dominated by it
//				if (strict == 1 && nd < dim) d = 1;

				else if (nd == dim) // dominates it
				{
					Erase.push_back(*it);
					outputFile3 << x[i] << ">=" << x[*it] << endl;
				}
			}

			if (!d) {
				for (auto it = Erase.begin(); it != Erase.end(); it++){ // see if we need to replace /excude any paretoindex dominated by i
					Paretoindex.erase(*it);
				}
				Paretoindex.insert(i);
//				cout <<endl<< i << " ";
//					for (int j = 0; j < 5; j++)  cout<<Data[dim*x[(i)] + j] << " " ;

			}


		}

		for (auto it = Paretoindex.begin(); it != Paretoindex.end(); it++) {
			d = 0;
			Erase.clear();

			for (auto it1 = Paretoindex.begin(); it1 != Paretoindex.end(); it1++) {
				nd = 0;
				if (*it != *it1) {
					if (strict == 1) { // strong pareto
						for (int j = 0;j < dim;j++) {
							if (Data[*it * dim + j] >= Data[*it1  *dim + j]) //not dominated    GB3 was >-
							{
								nd++;
							}
						}
					}
					else { // weak pareto
						for (int j = 0;j < dim;j++) {
							if (Data[*it *dim + j] > Data[*it1  *dim + j]) //not dominated    GB3 was >-
							{
								nd++;
							}
						}
					}
					if (nd == 0) d = 1; // dominated by it
	//				if (strict == 1 && nd < dim) d = 1;

					else if (nd == dim) // dominates it
					{
						Erase.push_back(*it1);
						outputFile3 << x[*it] << ">=" << x[*it1] << endl;
					}



				}
			}
				{
						for (auto it2 = Erase.begin(); it2 != Erase.end(); it2++) { // see if we need to replace /excude any paretoindex dominated by i
							Paretoindex.erase(*it2);
							outputFile3 << "erased " << x[*it2] <<  endl;
						}
						// already in Paretoindex.insert(i);
						//				cout <<endl<< i << " ";
						//					for (int j = 0; j < 5; j++)  cout<<Data[dim*x[(i)] + j] << " " ;

				}

		}
		outputFile1 << "#" << endl;
		for (auto it = Paretoindex.begin(); it != Paretoindex.end(); it++) {
			//outputFile1 << x[(*it )] << endl;
//			outputFile1 << *it << " ";
			for (int i = 0; i < 5; i++) outputFile1 << std::setprecision(12) <<Data[dim*(*it) + i]<<" ";   //GB2
			for (int i = 5; i < dim-1; i++) outputFile1 << std::setprecision(6) << Data[dim*(*it) + i] << " ";   //GB2v
			for (int i = dim-1; i < dim; i++) outputFile1 << std::setprecision(12) << Data[dim*(*it) + i] << " ";   //GB2v
			outputFile1 << endl;
			outputFile2 << x[(*it)] <<endl;
		}
		outputFile1 << "#" << endl;
		outputFile1.close();
		outputFile2.close();
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started:
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
