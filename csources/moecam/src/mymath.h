#ifndef MYMATH_H
#define MYMATH_H

#include <math.h>

// The TNT libary for array declarations and so on... 
#include "../tnt/tnt.h"

//using namespace TNT;					//breaking encapsulation of the TNT name space.
typedef TNT::Vector<int> iVec;			//type define iVec to hold integer vector 
typedef TNT::Vector<double> dVec;		//and dVec to hold double vectors.
typedef TNT::Matrix<double> dMatrix;	//holds matrixies of two or more dimentions.


const long double PI = 3.1415926535897932384626433832795;

//here are some functions written to handle some simple ops.
template<class NUMBER>
NUMBER MMdabs(NUMBER x)
{
	if(x>=0)
		return x;
	return x*(-1.0);
}

template<class NUMBER>
NUMBER MMmax(NUMBER x, NUMBER y)
{
	if(x>y)
		return x;
	return y;
}

template<class NUMBER>
NUMBER MMmin(NUMBER x, NUMBER y)
{
	if(x<y)
		return x;
	return y;
}

template<class DATATYPE>
int MMmin(TNT::Vector<DATATYPE>& x, int n)
{
	int i, i_min;
	DATATYPE min;

	for(i=2, i_min=1, min=x(1); i<=n; i++)
	{
		if(x(i)<min)
		{
			min=x(i);
			i_min=i;
		}
	}
	return i_min;
}

template<class DATATYPE>
int MMmax(TNT::Vector<DATATYPE>& x, int n)
{
	int i, i_max;
	DATATYPE max;

	for(i=2, i_max=1, max=x(1); i<=n; i++)
	{
		if(x(i)>max)
		{
			max=x(i);
			i_max=i;
		}
	}
	return i_max;
}

template<class DATATYPE>
int absmax(TNT::Vector<DATATYPE>& x, int n)
{
	int i, i_max;
	DATATYPE max;

	for(i=2, i_max=1, max=MMdabs(x(1)); i<=n; i++)
	{
		if(MMdabs(x(i))>max)
		{
			max=MMdabs(x(i));
			i_max=i;
		}
	}
	return i_max;
}

template<class DATATYPE>
void PrintVec(TNT::Vector<DATATYPE>& x, int n)
{
	int i;

	for(i=1; i<n; i++)
		cout<<x(i)<<",";
	cout<<x(i)<<endl;
}

template<class DATATYPE>
DATATYPE myRand(DATATYPE x=1, DATATYPE y=0)
{	 
	if(x>y)
		return y + ((x-y)*rand())/(RAND_MAX);
	else if(x<y)
		return x + ((y-x)*rand())/(RAND_MAX);

	return x;
}


//THIS ARE JUST SOME COMPARISON OPERATORS I MIGHT MAKE THEM MEMBERS OF THE CLUSTER CLASS
//OF COURSE THIS ONLY WORK FOR THOSE CLASSES WITH OVERLOADED OPERATOR== AND OPERATOR<
template<class DATATYPE>
bool operator< (const TNT::Vector<DATATYPE>& x, const TNT::Vector<DATATYPE>& y)
{
	int i, n;

	if(x.size() < y.size())
	{
		return false;
	}
	else if(x.size()==y.size())
	{
		n = x.size();
		for(i=1; i<=n; i++)
		{
			if(x(i)<y(i))
				return true;
			else if(x(i)==y(i))
				{;}
			else
				return false;
		}
	}

	return true;
}

template<class DATATYPE>
bool operator>= (const TNT::Vector<DATATYPE>& x, const TNT::Vector<DATATYPE>& y)
{
	return !(x<y);
}

template<class DATATYPE>
bool operator== (const TNT::Vector<DATATYPE>& x, const TNT::Vector<DATATYPE>& y)
{
	int i,n;

	if(x.size() != y.size())
	{
		return false;
	}
	
	n = x.size();
	for(i=1; i<=n; i++)
	{
		if(x(i)!=y(i))
			return false;
	}
	return true;
}

template<class DATATYPE>
bool operator!= (const TNT::Vector<DATATYPE>& x, const TNT::Vector<DATATYPE>& y)
{
	return !(x==y);
}

//THESE OPERATORS ARE REALLY SLOW MIGHT RE-WRITE IT.
template<class DATATYPE>
bool operator> (const TNT::Vector<DATATYPE>& x, const TNT::Vector<DATATYPE>& y)
{
	return !(x<y) && !(x==y);
}

template<class DATATYPE>
bool operator<= (const TNT::Vector<DATATYPE>& x, const TNT::Vector<DATATYPE>& y)
{
	return !(x>y);
}


//SOME OHTER FUNCTIONS I WILL NEED PLUSS SOME UNARY OPERATORS
template<class DATATYPE>
double sum(const TNT::Vector<DATATYPE>& x)
{
	int n = x.size(), i;
	DATATYPE sum =0;

	for(i=1;i<=n; i++)
		sum+=x(i);

	return sum;
}

template<class DATATYPE>
const TNT::Vector<DATATYPE> sqr(const TNT::Vector<DATATYPE>& x)
{
	int n = x.size(), i;
	TNT::Vector<DATATYPE> y;
	
	y.newsize(n);

	for(i=1; i<=n; i++)
		y(i)=x(i)*x(i);

	return y;
}

template<class DATATYPE>
const TNT::Vector<DATATYPE> sqrt(const TNT::Vector<DATATYPE>& x)
{
	int n = x.size(), i;
	TNT::Vector<DATATYPE> y;
	
	y.newsize(n);

	for(i=1; i<=n; i++)
		y(i)=sqrt(x(i));

	return y;
}

template<class DATATYPE>
const TNT::Vector<DATATYPE> operator/(const TNT::Vector<DATATYPE>& x, const double c)
{
	int n = x.size(), i;
	TNT::Vector<DATATYPE> y;
	
	y.newsize(n);

	for(i=1; i<=n; i++)
		y(i)=x(i)/c;

	return y;
}


#endif
