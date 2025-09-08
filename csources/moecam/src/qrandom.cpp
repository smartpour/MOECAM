#include "qrandom.h"

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

//# include "sobol.h"

//******************************************************************************

int i4_bit_hi1_base_2 ( int n )

//******************************************************************************
//
//  Purpose:
//
//    I4_BIT_HI1_BASE_2 returns the position of the high 1 bit base 2 in an integer.
//
//  Example:
//
//       N    Binary    Hi 1
//    ----    --------  ----
//       0           0     0
//       1           1     1
//       2          10     2
//       3          11     2 
//       4         100     3
//       5         101     3
//       6         110     3
//       7         111     3
//       8        1000     4
//       9        1001     4
//      10        1010     4
//      11        1011     4
//      12        1100     4
//      13        1101     4
//      14        1110     4
//      15        1111     4
//      16       10000     5
//      17       10001     5
//    1023  1111111111    10
//    1024 10000000000    11
//    1025 10000000001    11
//
//  Modified:
//
//    13 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the integer to be measured.
//    N should be nonnegative.  If N is nonpositive, I4_BIT_HI1_BASE_2
//    will always be 0.
//
//    Output, int I4_BIT_HI1_BASE_2, the number of bits base 2.
//
{
  int bit;

  bit = 0;

  while ( 0 < n )
  {
    bit = bit + 1;
    n = n / 2;
  }

  return bit;
}
//******************************************************************************

int i4_bit_lo0_base_2 ( int n )

//******************************************************************************
//
//  Purpose:
//
//    I4_BIT_LO0_BASE_2 returns the position of the low 0 bit base 2 in an integer.
//
//  Example:
//
//       N    Binary    Lo 0
//    ----    --------  ----
//       0           0     1
//       1           1     2
//       2          10     1
//       3          11     3 
//       4         100     1
//       5         101     2
//       6         110     1
//       7         111     4
//       8        1000     1
//       9        1001     2
//      10        1010     1
//      11        1011     3
//      12        1100     1
//      13        1101     2
//      14        1110     1
//      15        1111     5
//      16       10000     1
//      17       10001     2
//    1023  1111111111     1
//    1024 10000000000     1
//    1025 10000000001     1
//
//  Modified:
//
//    13 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the integer to be measured.
//    N should be nonnegative.
//
//    Output, int I4_BIT_LO0_BASE_2, the position of the low 1 bit.
//
{
  int bit;
  int n2;

  bit = 0;

  while ( true )
  {
    bit = bit + 1;
    n2 = n / 2;

    if ( n == 2 * n2 )
    {
      break;
    }

    n = n2;

  }

  return bit;
}
//******************************************************************************

void i4_sobol ( int dim_num, int *seed, float quasi[ ] )

//******************************************************************************
//
//  Purpose:
//
//    I4_SOBOL generates a new quasirandom Sobol vector with each call.
//
//  Discussion:
//
//    The routine adapts the ideas of Antonov and Saleev.
//
//  Modified:
//
//    03 August 2004
//
//  Reference:
//
//    Antonov and Saleev,
//    USSR Computational Mathematics and Mathematical Physics,
//    Volume 19, 1980, pages 252 - 256.
//
//    Paul Bratley and Bennett Fox,
//    Algorithm 659:
//    Implementing Sobol's Quasirandom Sequence Generator,
//    ACM Transactions on Mathematical Software,
//    Volume 14, Number 1, pages 88-100, 1988.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom 
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    I Sobol,
//    USSR Computational Mathematics and Mathematical Physics,
//    Volume 16, pages 236-242, 1977.
//
//    I Sobol and Levitan, 
//    The Production of Points Uniformly Distributed in a Multidimensional 
//    Cube (in Russian),
//    Preprint IPM Akad. Nauk SSSR, 
//    Number 40, Moscow 1976.
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//    DIM_NUM must satisfy 2 <= DIM_NUM <= 40.
//
//    Input/output, int *SEED, the "seed" for the sequence.
//    This is essentially the index in the sequence of the quasirandom
//    value to be generated.  On output, SEED has been set to the
//    appropriate next value, usually simply SEED+1.
//    If SEED is less than 0 on input, it is treated as though it were 0.
//    An input value of 0 requests the first (0-th) element of the sequence.
//
//    Output, float QUASI(DIM_NUM), the next quasirandom vector.
//
{
# define DIM_MAX 40

  static int atmost = 1073741823;
  static int dim_num_save = 0;
  int i;
  int i2;
  bool includ[8];
  static bool initialized = false;
  int j;
  int j2;
  int k;
  int l;
  static int lastq[DIM_MAX];
  int m;
  static int maxcol;
  int newv;
  static int poly[DIM_MAX] =
  {
        1,   3,   7,  11,  13,  19,  25,  37,  59,  47,
       61,  55,  41,  67,  97,  91, 109, 103, 115, 131,
      193, 137, 145, 143, 241, 157, 185, 167, 229, 171,
      213, 191, 253, 203, 211, 239, 247, 285, 369, 299 
  };
  static float recipd;
  static int seed_save = 0;
  int seed_temp;
  static int v[DIM_MAX][30];
//
  if ( !initialized || dim_num != dim_num_save )
  {
    initialized = true;
//
//  Initialize (part of) V.
//
    v[ 0][0] = 1;
    v[ 1][0] = 1;
    v[ 2][0] = 1;
    v[ 3][0] = 1;
    v[ 4][0] = 1;
    v[ 5][0] = 1;
    v[ 6][0] = 1;
    v[ 7][0] = 1;
    v[ 8][0] = 1;
    v[ 9][0] = 1;
    v[10][0] = 1;
    v[11][0] = 1;
    v[12][0] = 1;
    v[13][0] = 1;
    v[14][0] = 1;
    v[15][0] = 1;
    v[16][0] = 1;
    v[17][0] = 1;
    v[18][0] = 1;
    v[19][0] = 1;
    v[20][0] = 1;
    v[21][0] = 1;
    v[22][0] = 1;
    v[23][0] = 1;
    v[24][0] = 1;
    v[25][0] = 1;
    v[26][0] = 1;
    v[27][0] = 1;
    v[28][0] = 1;
    v[29][0] = 1;
    v[30][0] = 1;
    v[31][0] = 1;
    v[32][0] = 1;
    v[33][0] = 1;
    v[34][0] = 1;
    v[35][0] = 1;
    v[36][0] = 1;
    v[37][0] = 1;
    v[38][0] = 1;
    v[39][0] = 1;

    v[ 2][1] = 1;
    v[ 3][1] = 3;
    v[ 4][1] = 1;
    v[ 5][1] = 3;
    v[ 6][1] = 1;
    v[ 7][1] = 3;
    v[ 8][1] = 3;
    v[ 9][1] = 1;
    v[10][1] = 3;
    v[11][1] = 1;
    v[12][1] = 3;
    v[13][1] = 1;
    v[14][1] = 3;
    v[15][1] = 1;
    v[16][1] = 1;
    v[17][1] = 3;
    v[18][1] = 1;
    v[19][1] = 3;
    v[20][1] = 1;
    v[21][1] = 3;
    v[22][1] = 1;
    v[23][1] = 3;
    v[24][1] = 3;
    v[25][1] = 1;
    v[26][1] = 3;
    v[27][1] = 1;
    v[28][1] = 3;
    v[29][1] = 1;
    v[30][1] = 3;
    v[31][1] = 1;
    v[32][1] = 1;
    v[33][1] = 3;
    v[34][1] = 1;
    v[35][1] = 3;
    v[36][1] = 1;
    v[37][1] = 3;
    v[38][1] = 1;
    v[39][1] = 3;

    v[ 3][2] = 7;
    v[ 4][2] = 5;
    v[ 5][2] = 1;
    v[ 6][2] = 3;
    v[ 7][2] = 3;
    v[ 8][2] = 7;
    v[ 9][2] = 5;
    v[10][2] = 5;
    v[11][2] = 7;
    v[12][2] = 7;
    v[13][2] = 1;
    v[14][2] = 3;
    v[15][2] = 3;
    v[16][2] = 7;
    v[17][2] = 5;
    v[18][2] = 1;
    v[19][2] = 1;
    v[20][2] = 5;
    v[21][2] = 3;
    v[22][2] = 3;
    v[23][2] = 1;
    v[24][2] = 7;
    v[25][2] = 5;
    v[26][2] = 1;
    v[27][2] = 3;
    v[28][2] = 3;
    v[29][2] = 7;
    v[30][2] = 5;
    v[31][2] = 1;
    v[32][2] = 1;
    v[33][2] = 5;
    v[34][2] = 7;
    v[35][2] = 7;
    v[36][2] = 5;
    v[37][2] = 1;
    v[38][2] = 3;
    v[39][2] = 3;

    v[ 5][3] =  1;
    v[ 6][3] =  7;
    v[ 7][3] =  9;
    v[ 8][3] = 13;
    v[ 9][3] = 11;
    v[10][3] =  1;
    v[11][3] =  3;
    v[12][3] =  7;
    v[13][3] =  9;
    v[14][3] =  5;
    v[15][3] = 13;
    v[16][3] = 13;
    v[17][3] = 11;
    v[18][3] =  3;
    v[19][3] = 15;
    v[20][3] =  5;
    v[21][3] =  3;
    v[22][3] = 15;
    v[23][3] =  7;
    v[24][3] =  9;
    v[25][3] = 13;
    v[26][3] =  9;
    v[27][3] =  1;
    v[28][3] = 11;
    v[29][3] =  7;
    v[30][3] =  5;
    v[31][3] = 15;
    v[32][3] =  1;
    v[33][3] = 15;
    v[34][3] = 11;
    v[35][3] =  5;
    v[36][3] =  3;
    v[37][3] =  1;
    v[38][3] =  7;
    v[39][3] =  9;
  
    v[ 7][4] =  9;
    v[ 8][4] =  3;
    v[ 9][4] = 27;
    v[10][4] = 15;
    v[11][4] = 29;
    v[12][4] = 21;
    v[13][4] = 23;
    v[14][4] = 19;
    v[15][4] = 11;
    v[16][4] = 25;
    v[17][4] =  7;
    v[18][4] = 13;
    v[19][4] = 17;
    v[20][4] =  1;
    v[21][4] = 25;
    v[22][4] = 29;
    v[23][4] =  3;
    v[24][4] = 31;
    v[25][4] = 11;
    v[26][4] =  5;
    v[27][4] = 23;
    v[28][4] = 27;
    v[29][4] = 19;
    v[30][4] = 21;
    v[31][4] =  5;
    v[32][4] =  1;
    v[33][4] = 17;
    v[34][4] = 13;
    v[35][4] =  7;
    v[36][4] = 15;
    v[37][4] =  9;
    v[38][4] = 31;
    v[39][4] =  9;

    v[13][5] = 37;
    v[14][5] = 33;
    v[15][5] =  7;
    v[16][5] =  5;
    v[17][5] = 11;
    v[18][5] = 39;
    v[19][5] = 63;
    v[20][5] = 27;
    v[21][5] = 17;
    v[22][5] = 15;
    v[23][5] = 23;
    v[24][5] = 29;
    v[25][5] =  3;
    v[26][5] = 21;
    v[27][5] = 13;
    v[28][5] = 31;
    v[29][5] = 25;
    v[30][5] =  9;
    v[31][5] = 49;
    v[32][5] = 33;
    v[33][5] = 19;
    v[34][5] = 29;
    v[35][5] = 11;
    v[36][5] = 19;
    v[37][5] = 27;
    v[38][5] = 15;
    v[39][5] = 25;

    v[19][6] =  13;
    v[20][6] =  35;
    v[21][6] = 115;
    v[22][6] =  41;
    v[23][6] =  79;
    v[24][6] =  17;
    v[25][6] =  29;
    v[26][6] = 119;
    v[27][6] =  75;
    v[28][6] =  73;
    v[29][6] = 105;
    v[30][6] =   7;
    v[31][6] =  59;
    v[32][6] =  65;
    v[33][6] =  21;
    v[34][6] =   3;
    v[35][6] = 113;
    v[36][6] =  61;
    v[37][6] =  89;
    v[38][6] =  45;
    v[39][6] = 107;

    v[37][7] =  7;
    v[38][7] = 23;
    v[39][7] = 39;
//
//  Check parameters.
//
    if ( dim_num < 2 || DIM_MAX < dim_num )
    {
      cout << "\n";
      cout << "I4_SOBOL - Fatal error!\n";
      cout << "  The spatial dimension DIM_NUM should satisfy:\n";
      cout << "    2 <= DIM_NUM <= " << DIM_MAX << "\n";
      cout << "  But this input value is DIM_NUM = " << dim_num << "\n";
      exit ( 1 );
    }

    dim_num_save = dim_num;
//
//  Find the number of bits in ATMOST.
//
    maxcol = i4_bit_hi1_base_2 ( atmost );
//
//  Initialize row 1 of V.
//
    for ( j = 1; j <= maxcol; j++ )
    {
      v[1-1][j-1] = 1;
    }
//
//  Initialize the remaining rows of V.
//
    for ( i = 1; i < dim_num; i++ )
    {
//
//  The bit pattern of the integer POLY(I) gives the form
//  of polynomial I.
//
//  Find the degree of polynomial I from binary encoding.
//
      j = poly[i];
      m = 0;

      while ( true )
      {
        j = j / 2;
        if ( j <= 0 )
        {
          break;
        }
        m = m + 1;
      }
//
//  We expand this bit pattern to separate components
//  of the logical array INCLUD.
//
      j = poly[i];
      for ( k = m-1; k >= 0; k-- )
      {
        j2 = j / 2;
        includ[k] = ( j != ( 2 * j2 ) );
        j = j2;
      }
//
//  Calculate the remaining elements of row I as explained
//  in Bratley and Fox, section 2.
//
//  Some tricky indexing here.  Did I change it correctly?
//
      for ( j = m; j < maxcol; j++ )
      {
        newv = v[i][j-m];
        l = 1;

        for ( k = 0; k < m; k++ )
        {
          l = 2 * l;

          if ( includ[k] )
          {
            newv = ( newv ^ ( l * v[i][j-k-1] ) );
          }

        }

        v[i][j] = newv;

      }

    }
//
//  Multiply columns of V by appropriate power of 2.
//
    l = 1;
    for ( j = maxcol-2; j >= 0; j-- )
    {
      l = 2 * l;
      for ( i = 0; i < dim_num; i++ )
      {
        v[i][j] = v[i][j] * l;
      }
    }
//
//  RECIPD is 1/(common denominator of the elements in V).
//
    recipd = 1.0E+00 / ( ( float ) ( 2 * l ) );
  }

  if ( *seed < 0 )
  {
    *seed = 0;
  }

  if ( *seed == 0 )
  {
    l = 1;
    for ( i = 0; i < dim_num; i++ )
    {
      lastq[i] = 0;
    }
  }
  else if ( *seed == seed_save + 1 )
  {
    l = i4_bit_lo0_base_2 ( *seed );
  }
  else if ( *seed <= seed_save )
  {
    seed_save = 0;
    l = 1;
    for ( i = 0; i < dim_num; i++ )
    {
      lastq[i] = 0;
    }

    for ( seed_temp = seed_save; seed_temp <= (*seed)-1; seed_temp++ )
    {

      l = i4_bit_lo0_base_2 ( seed_temp );

      for ( i = 0; i < dim_num; i++ )
      {
        lastq[i] = ( lastq[i] ^ v[i][l-1] );
      }

    }

    l = i4_bit_lo0_base_2 ( *seed );
  }
  else if ( seed_save+1 < *seed )
  {
    for ( seed_temp = seed_save+1; seed_temp <= (*seed)-1; seed_temp++ )
    {

      l = i4_bit_lo0_base_2 ( seed_temp );

      for ( i = 0; i < dim_num; i++ )
      {
        lastq[i] = ( lastq[i] ^ v[i][l-1] );
      }

    }

    l = i4_bit_lo0_base_2 ( *seed );

  }
//
//  Check that the user is not calling too many times!
//
  if ( maxcol < l )
  {
    cout << "\n";
    cout << "I4_SOBOL - Fatal error!\n";
    cout << "  Too many calls!\n";
    cout << "  MAXCOL = " << maxcol << "\n";
    cout << "  L =      " << l << "\n";
    exit ( 2 );
  }
//
//  Calculate the new components of QUASI.
//  The caret indicates the bitwise exclusive OR.
//
  for ( i = 0; i < dim_num; i++ )
  {
    quasi[i] = ( ( float ) lastq[i] ) * recipd;

    lastq[i] = ( lastq[i] ^ v[i][l-1] );
  }

  seed_save = *seed;
  *seed = *seed + 1;

  return;
# undef MAX_DIM
}
//******************************************************************************

int i8_bit_hi1_base_2 (unsigned long int n )

//******************************************************************************
//
//  Purpose:
//
//    I8_BIT_HI1_BASE_2 returns the position of the high 1 bit base 2 in an integer.
//
//  Example:
//
//       N    Binary    Hi 1
//    ----    --------  ----
//       0           0     0
//       1           1     1
//       2          10     2
//       3          11     2 
//       4         100     3
//       5         101     3
//       6         110     3
//       7         111     3
//       8        1000     4
//       9        1001     4
//      10        1010     4
//      11        1011     4
//      12        1100     4
//      13        1101     4
//      14        1110     4
//      15        1111     4
//      16       10000     5
//      17       10001     5
//    1023  1111111111    10
//    1024 10000000000    11
//    1025 10000000001    11
//
//  Modified:
//
//    03 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, long int N, the integer to be measured.
//    N should be nonnegative.  If N is nonpositive, I8_BIT_HI1_BASE_2
//    will always be 0.
//
//    Output, int I8_BIT_HI1_BASE_2, the number of bits base 2.
//
{
  int bit;

  bit = 0;

  while ( 0 < n )
  {
    bit = bit + 1;
    n = n / 2;
  }

  return bit;
}
//******************************************************************************

int i8_bit_lo0_base_2 ( long int n )

//******************************************************************************
//
//  Purpose:
//
//    I8_BIT_LO0_BASE_2 returns the position of the low 0 bit base 2 in an integer.
//
//  Example:
//
//       N    Binary    Lo 0
//    ----    --------  ----
//       0           0     1
//       1           1     2
//       2          10     1
//       3          11     3 
//       4         100     1
//       5         101     2
//       6         110     1
//       7         111     4
//       8        1000     1
//       9        1001     2
//      10        1010     1
//      11        1011     3
//      12        1100     1
//      13        1101     2
//      14        1110     1
//      15        1111     5
//      16       10000     1
//      17       10001     2
//    1023  1111111111     1
//    1024 10000000000     1
//    1025 10000000001     1
//
//  Modified:
//
//    03 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, long int N, the integer to be measured.
//    N should be nonnegative.
//
//    Output, int I8_BIT_LO0_BASE_2, the position of the low 1 bit.
//
{
  int bit;
  long int n2;

  bit = 0;

  while ( true )
  {
    bit = bit + 1;
    n2 = n / 2;

    if ( n == 2 * n2 )
    {
      break;
    }

    n = n2;

  }

  return bit;
}
//******************************************************************************

void i8_sobol ( int dim_num, long int *seed, double quasi[ ] )

//******************************************************************************
//
//  Purpose:
//
//    I8_SOBOL generates a new quasirandom Sobol vector with each call.
//
//  Discussion:
//
//    The routine adapts the ideas of Antonov and Saleev.
//
//  Modified:
//
//    03 August 2004
//
//  Reference:
//
//    Antonov and Saleev,
//    USSR Computational Mathematics and Mathematical Physics,
//    Volume 19, 1980, pages 252 - 256.
//
//    Paul Bratley and Bennett Fox,
//    Algorithm 659:
//    Implementing Sobol's Quasirandom Sequence Generator,
//    ACM Transactions on Mathematical Software,
//    Volume 14, Number 1, pages 88-100, 1988.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom 
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    I Sobol,
//    USSR Computational Mathematics and Mathematical Physics,
//    Volume 16, pages 236-242, 1977.
//
//    I Sobol and Levitan, 
//    The Production of Points Uniformly Distributed in a Multidimensional 
//    Cube (in Russian),
//    Preprint IPM Akad. Nauk SSSR, 
//    Number 40, Moscow 1976.
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//    DIM_NUM must satisfy 2 <= DIM_NUM <= 40.
//
//    Input/output, long int *SEED, the "seed" for the sequence.
//    This is essentially the index in the sequence of the quasirandom
//    value to be generated.  On output, SEED has been set to the
//    appropriate next value, usually simply SEED+1.
//    If SEED is less than 0 on input, it is treated as though it were 0.
//    An input value of 0 requests the first (0-th) element of the sequence.
//
//    Output, double QUASI[DIM_NUM], the next quasirandom vector.
//
{
# define DIM_MAX 40

  //static unsigned long int atmost = 4611686018427387903ULL;
  static unsigned long int atmost = 1073741823UL;
  static int dim_num_save = 0;
  int i;
  int i2;
  bool includ[8];
  static bool initialized = false;
  int j;
  int j2;
  int k;
  long int l;
  static long int lastq[DIM_MAX];
  int m;
  static int maxcol;
  long int newv;
  static int poly[DIM_MAX] =
  {
        1,   3,   7,  11,  13,  19,  25,  37,  59,  47,
       61,  55,  41,  67,  97,  91, 109, 103, 115, 131,
      193, 137, 145, 143, 241, 157, 185, 167, 229, 171,
      213, 191, 253, 203, 211, 239, 247, 285, 369, 299 
  };
  static double recipd;
  static long int seed_save = 0;
  long int seed_temp;
  static long int v[DIM_MAX][30];
//
  if ( !initialized || dim_num != dim_num_save )
  {
    initialized = true;
//
//  Initialize (part of) V.
//
    v[ 0][0] = 1;
    v[ 1][0] = 1;
    v[ 2][0] = 1;
    v[ 3][0] = 1;
    v[ 4][0] = 1;
    v[ 5][0] = 1;
    v[ 6][0] = 1;
    v[ 7][0] = 1;
    v[ 8][0] = 1;
    v[ 9][0] = 1;
    v[10][0] = 1;
    v[11][0] = 1;
    v[12][0] = 1;
    v[13][0] = 1;
    v[14][0] = 1;
    v[15][0] = 1;
    v[16][0] = 1;
    v[17][0] = 1;
    v[18][0] = 1;
    v[19][0] = 1;
    v[20][0] = 1;
    v[21][0] = 1;
    v[22][0] = 1;
    v[23][0] = 1;
    v[24][0] = 1;
    v[25][0] = 1;
    v[26][0] = 1;
    v[27][0] = 1;
    v[28][0] = 1;
    v[29][0] = 1;
    v[30][0] = 1;
    v[31][0] = 1;
    v[32][0] = 1;
    v[33][0] = 1;
    v[34][0] = 1;
    v[35][0] = 1;
    v[36][0] = 1;
    v[37][0] = 1;
    v[38][0] = 1;
    v[39][0] = 1;

    v[ 2][1] = 1;
    v[ 3][1] = 3;
    v[ 4][1] = 1;
    v[ 5][1] = 3;
    v[ 6][1] = 1;
    v[ 7][1] = 3;
    v[ 8][1] = 3;
    v[ 9][1] = 1;
    v[10][1] = 3;
    v[11][1] = 1;
    v[12][1] = 3;
    v[13][1] = 1;
    v[14][1] = 3;
    v[15][1] = 1;
    v[16][1] = 1;
    v[17][1] = 3;
    v[18][1] = 1;
    v[19][1] = 3;
    v[20][1] = 1;
    v[21][1] = 3;
    v[22][1] = 1;
    v[23][1] = 3;
    v[24][1] = 3;
    v[25][1] = 1;
    v[26][1] = 3;
    v[27][1] = 1;
    v[28][1] = 3;
    v[29][1] = 1;
    v[30][1] = 3;
    v[31][1] = 1;
    v[32][1] = 1;
    v[33][1] = 3;
    v[34][1] = 1;
    v[35][1] = 3;
    v[36][1] = 1;
    v[37][1] = 3;
    v[38][1] = 1;
    v[39][1] = 3;

    v[ 3][2] = 7;
    v[ 4][2] = 5;
    v[ 5][2] = 1;
    v[ 6][2] = 3;
    v[ 7][2] = 3;
    v[ 8][2] = 7;
    v[ 9][2] = 5;
    v[10][2] = 5;
    v[11][2] = 7;
    v[12][2] = 7;
    v[13][2] = 1;
    v[14][2] = 3;
    v[15][2] = 3;
    v[16][2] = 7;
    v[17][2] = 5;
    v[18][2] = 1;
    v[19][2] = 1;
    v[20][2] = 5;
    v[21][2] = 3;
    v[22][2] = 3;
    v[23][2] = 1;
    v[24][2] = 7;
    v[25][2] = 5;
    v[26][2] = 1;
    v[27][2] = 3;
    v[28][2] = 3;
    v[29][2] = 7;
    v[30][2] = 5;
    v[31][2] = 1;
    v[32][2] = 1;
    v[33][2] = 5;
    v[34][2] = 7;
    v[35][2] = 7;
    v[36][2] = 5;
    v[37][2] = 1;
    v[38][2] = 3;
    v[39][2] = 3;

    v[ 5][3] =  1;
    v[ 6][3] =  7;
    v[ 7][3] =  9;
    v[ 8][3] = 13;
    v[ 9][3] = 11;
    v[10][3] =  1;
    v[11][3] =  3;
    v[12][3] =  7;
    v[13][3] =  9;
    v[14][3] =  5;
    v[15][3] = 13;
    v[16][3] = 13;
    v[17][3] = 11;
    v[18][3] =  3;
    v[19][3] = 15;
    v[20][3] =  5;
    v[21][3] =  3;
    v[22][3] = 15;
    v[23][3] =  7;
    v[24][3] =  9;
    v[25][3] = 13;
    v[26][3] =  9;
    v[27][3] =  1;
    v[28][3] = 11;
    v[29][3] =  7;
    v[30][3] =  5;
    v[31][3] = 15;
    v[32][3] =  1;
    v[33][3] = 15;
    v[34][3] = 11;
    v[35][3] =  5;
    v[36][3] =  3;
    v[37][3] =  1;
    v[38][3] =  7;
    v[39][3] =  9;
  
    v[ 7][4] =  9;
    v[ 8][4] =  3;
    v[ 9][4] = 27;
    v[10][4] = 15;
    v[11][4] = 29;
    v[12][4] = 21;
    v[13][4] = 23;
    v[14][4] = 19;
    v[15][4] = 11;
    v[16][4] = 25;
    v[17][4] =  7;
    v[18][4] = 13;
    v[19][4] = 17;
    v[20][4] =  1;
    v[21][4] = 25;
    v[22][4] = 29;
    v[23][4] =  3;
    v[24][4] = 31;
    v[25][4] = 11;
    v[26][4] =  5;
    v[27][4] = 23;
    v[28][4] = 27;
    v[29][4] = 19;
    v[30][4] = 21;
    v[31][4] =  5;
    v[32][4] =  1;
    v[33][4] = 17;
    v[34][4] = 13;
    v[35][4] =  7;
    v[36][4] = 15;
    v[37][4] =  9;
    v[38][4] = 31;
    v[39][4] =  9;

    v[13][5] = 37;
    v[14][5] = 33;
    v[15][5] =  7;
    v[16][5] =  5;
    v[17][5] = 11;
    v[18][5] = 39;
    v[19][5] = 63;
    v[20][5] = 27;
    v[21][5] = 17;
    v[22][5] = 15;
    v[23][5] = 23;
    v[24][5] = 29;
    v[25][5] =  3;
    v[26][5] = 21;
    v[27][5] = 13;
    v[28][5] = 31;
    v[29][5] = 25;
    v[30][5] =  9;
    v[31][5] = 49;
    v[32][5] = 33;
    v[33][5] = 19;
    v[34][5] = 29;
    v[35][5] = 11;
    v[36][5] = 19;
    v[37][5] = 27;
    v[38][5] = 15;
    v[39][5] = 25;

    v[19][6] =  13;
    v[20][6] =  35;
    v[21][6] = 115;
    v[22][6] =  41;
    v[23][6] =  79;
    v[24][6] =  17;
    v[25][6] =  29;
    v[26][6] = 119;
    v[27][6] =  75;
    v[28][6] =  73;
    v[29][6] = 105;
    v[30][6] =   7;
    v[31][6] =  59;
    v[32][6] =  65;
    v[33][6] =  21;
    v[34][6] =   3;
    v[35][6] = 113;
    v[36][6] =  61;
    v[37][6] =  89;
    v[38][6] =  45;
    v[39][6] = 107;

    v[37][7] =  7;
    v[38][7] = 23;
    v[39][7] = 39;
//
//  Check parameters.
//
    if ( dim_num < 2 || DIM_MAX < dim_num )
    {
      cout << "\n";
      cout << "I8_SOBOL - Fatal error!\n";
      cout << "  The spatial dimension DIM_NUM should satisfy:\n";
      cout << "    2 <= DIM_NUM <= " << DIM_MAX << "\n";
      cout << "  But this input value is DIM_NUM = " << dim_num << "\n";
      exit ( 1 );
    }

    dim_num_save = dim_num;
//
//  Find the number of bits in ATMOST.
//
    maxcol = i8_bit_hi1_base_2 ( atmost );
//
//  Initialize row 1 of V.
//
    for ( j = 1; j <= maxcol; j++ )
    {
      v[1-1][j-1] = 1;
    }
//
//  Initialize the remaining rows of V.
//
    for ( i = 1; i < dim_num; i++ )
    {
//
//  The bit pattern of the integer POLY(I) gives the form
//  of polynomial I.
//
//  Find the degree of polynomial I from binary encoding.
//
      j = poly[i];
      m = 0;

      while ( true )
      {
        j = j / 2;
        if ( j <= 0 )
        {
          break;
        }
        m = m + 1;
      }
//
//  We expand this bit pattern to separate components
//  of the logical array INCLUD.
//
      j = poly[i];
      for ( k = m-1; k >= 0; k-- )
      {
        j2 = j / 2;
        includ[k] = ( j != ( 2 * j2 ) );
        j = j2;
      }
//
//  Calculate the remaining elements of row I as explained
//  in Bratley and Fox, section 2.
//
//  Some tricky indexing here.  Did I change it correctly?
//
      for ( j = m; j < maxcol; j++ )
      {
        newv = v[i][j-m];
        l = 1;

        for ( k = 0; k < m; k++ )
        {
          l = 2 * l;

          if ( includ[k] )
          {
            newv = ( newv ^ ( l * v[i][j-k-1] ) );
          }

        }

        v[i][j] = newv;

      }

    }
//
//  Multiply columns of V by appropriate power of 2.
//
    l = 1;
    for ( j = maxcol-2; j >= 0; j-- )
    {
      l = 2 * l;
      for ( i = 0; i < dim_num; i++ )
      {
        v[i][j] = v[i][j] * l;
      }
    }
//
//  RECIPD is 1/(common denominator of the elements in V).
//
    recipd = 1.0E+00 / ( ( double ) ( 2 * l ) );
  }

  if ( *seed < 0 )
  {
    *seed = 0;
  }

  if ( *seed == 0 )
  {
    l = 1;
    for ( i = 0; i < dim_num; i++ )
    {
      lastq[i] = 0;
    }
  }
  else if ( *seed == seed_save + 1 )
  {
    l = i8_bit_lo0_base_2 ( *seed );
  }
  else if ( *seed <= seed_save )
  {
    seed_save = 0;
    l = 1;
    for ( i = 0; i < dim_num; i++ )
    {
      lastq[i] = 0;
    }

    for ( seed_temp = seed_save; seed_temp <= (*seed)-1; seed_temp++ )
    {

      l = i8_bit_lo0_base_2 ( seed_temp );

      for ( i = 0; i < dim_num; i++ )
      {
        lastq[i] = ( lastq[i] ^ v[i][l-1] );
      }

    }

    l = i8_bit_lo0_base_2 ( *seed );
  }
  else if ( seed_save+1 < *seed )
  {
    for ( seed_temp = seed_save+1; seed_temp <= (*seed)-1; seed_temp++ )
    {

      l = i8_bit_lo0_base_2 ( seed_temp );

      for ( i = 0; i < dim_num; i++ )
      {
        lastq[i] = ( lastq[i] ^ v[i][l-1] );
      }

    }

    l = i8_bit_lo0_base_2 ( *seed );

  }
//
//  Check that the user is not calling too many times!
//
  if ( maxcol < l )
  {
    cout << "\n";
    cout << "I8_SOBOL - Fatal error!\n";
    cout << "  Too many calls!\n";
    cout << "  MAXCOL = " << maxcol << "\n";
    cout << "  L =      " << l << "\n";
    exit ( 2 );
  }
//
//  Calculate the new components of QUASI.
//  The caret indicates the bitwise exclusive OR.
//
  for ( i = 0; i < dim_num; i++ )
  {
    quasi[i] = ( ( double ) lastq[i] ) * recipd;

    lastq[i] = ( lastq[i] ^ v[i][l-1] );
  }

  seed_save = *seed;
  *seed = *seed + 1;

  return;
# undef MAX_DIM
}
//**********************************************************************

void timestamp ( void )

//**********************************************************************
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Modified:
//
//    04 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
#define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
#undef TIME_SIZE
}


#include <vector>
#include <algorithm>
void	qrandom::RandomVec(int dim, double* x)
{
	if(dim<40)
		i8_sobol ( dim, &seed, x );
	else 
		RandomUniform(dim,x);
}

void	qrandom::RandomUniform(int dim, double* x)
{
	int i;
	for(i=0;i<dim;i++) x[i]=rand()/double(RAND_MAX);
}


void	qrandom::RandomVecSimplex(int dim, double * x)
{
	int i;
	if((dim-1)<40)
		i8_sobol ( dim-1, &seed, x );
	else RandomUniform(dim-1,x);
	x[dim-1]=1;
	std::sort(&x[0],&x[dim-1]); // increasing order
	for(i=dim-1;i>0;i--) x[i]=x[i]-x[i-1];
}


