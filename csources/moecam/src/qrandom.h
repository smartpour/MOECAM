/*
	interface to quasirandom sequence generator (Sobol sequence)
*/



#define MAX_DIM 40
class qrandom {
public:
	long int seed;
	qrandom() {seed=1;}

	// uniform quasirandom vector of size 1 < dim <= Max_DIM
	void	RandomVec(int dim, double* x);

	// uniform quasirandom vector of size 1 < dim <= Max_DIM on the unit simplex
	// sum x[i] =1  dim is the dimension of the simplex (ie 1+dim of the space it is in)
	void	RandomVecSimplex(int dim, double * x);

	// standard C rand() call for dim>40
	void	RandomUniform(int dim, double* x);


};
