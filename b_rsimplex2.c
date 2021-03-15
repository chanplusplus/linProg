#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "vector.h"


double abs_max( double a, double b )
{
	if( fabs( a ) > fabs( b ) )
		return a;
	else
		return b;
}

void solution( int m, int n, double **Binv, int *basis, double *x )
{
	int i;
	for( i = 0; i < n; i++ ) x[i] = 0;
	for( i = 0; i < m; i++ )
	{
		if( basis[i] < n ) x[basis[i]] = Binv[i][m];
	}
}


void simplex_multiplier( int m, double **Binv, double *b, double *c, int *basis )
{
	int i,j;
	for( i = 0; i < m; i++ )
	{
		Binv[m][i] = 0;
		for( j = 0; j < m; j++ )
		{
			Binv[m][i] = Binv[m][i] - c[basis[j]] * Binv[j][i];
		}
	}
	Binv[m][m] = inner( Binv[m], b, m );
}


void pivoting( int m, double **Binv, double *p_hat, double c_s, int r )
{
	int i, j;
	for( i = 0; i < m; i++ )
	{
		if( i != r )
			for( j = 0; j < m + 1; j++ ) Binv[i][j] = Binv[i][j] - Binv[r][j] * p_hat[i] / p_hat[r];
	}
	for( j = 0; j < m + 1; j++ ) Binv[m][j] = Binv[m][j] - Binv[r][j] * c_s / p_hat[r];
	for( j = 0; j < m + 1; j++ ) Binv[r][j] = Binv[r][j] / p_hat[r];
}

void tableau( int m, int n, double **Binv, double **p, double **A, double *b, double *c, double *d, int *basis, int *basisset )
{
	int i, j;



	for( i = 0; i < m; i++ )
		for( j = 0; j < n; j++ ) p[j][i] = A[i][j];
	for( j = 0; j < n; j++ ) p[j][m] = c[j];



	for( i = 0; i < m; i++ ) Binv[i][m+1] = b[i] * ( double )( ( b[i] > -EPS ) ? 1 : -1 );
	for( i = 0; i < m + 1; i++ ) Binv[i][i] = 1;
	for( i = 0; i < m; i++ ) p[n+i][i] = 1;
	for( i = 0; i < m; i++ )
	{
		if( b[i] > -EPS )
		{
			basis[i] = n + i;
			basisset[n+i] = 1;
		}
		else{
			basis[i] = n + m;
			for( j = 0; j < n + m; j++ )
			{
				p[j][i] = -p[j][i];
				d[j] = d[j] - p[j][i];
			}
			Binv[m+1][m+1] = Binv[m+1][m+1] - Binv[i][m+1];
		}
	}
}


void find_s( int m, int n, double **Binv, double **p, double *c, int *basisset, int *s, double *c_s )
{
	int j;
	double temp;
	*c_s = INFTY;
  
	for( j = 0; j < n; j++ )
	{
		if( basisset[j] == 0 )
		{
			temp = c[j] + inner( Binv[m], p[j], m );
			if( *c_s > temp )
			{
				*c_s = temp;
				*s = j;
			}
		}
	}
}



void b_find_s( int m, int n, double **Binv, double **p, double *c, int *basisset, int *s, double *c_s )
{
	int j;
	double temp;
	*c_s = INFTY;

	for( j = 0; j < n; j++ )
	{
		if( basisset[j] == 0 )
		{
			temp = c[j] + inner( Binv[m], p[j], m );
			if( temp < -EPS )
			{
				*c_s = temp;
				*s = j;
				return;
			}
		}
	}
}


void r_find_s( int m, int n, double **Binv, double **p, double *c, double *p_hat, int *basisset, int *s, double *c_s )
{
	int i, j;
	double temp;
	double temp_ratio = INFTY;
	*c_s = INFTY;

	for( j = 0; j < n; j++ )
	{
		if( basisset[j] == 0 )
		{
			for( i = 0; i < m; i++ ) p_hat[i] = inner( Binv[i], p[j], m );

			for( i = 0; i < m; i++ )
			{
				if( p_hat[i] > EPS && Binv[i][m] > -EPS )
				{
					if( temp_ratio > Binv[i][m] / p_hat[i] )
					{
						temp_ratio = Binv[i][m] / p_hat[i];
					}
				}
			}

			temp = ( c[j] + inner( Binv[m], p[j], m ) ) * temp_ratio;
			if( *c_s > temp )
			{
				*c_s = temp;
				*s = j;
			}
		}
	}
}


void find_r( int m, int n, double **Binv, double *p_hat, int *r )
{
	int i;
	double temp;

	temp = INFTY;
	for( i = 0; i < n; i++ )
	{
		if( p_hat[i] > EPS && Binv[i][m] > -EPS )
		{
			if( temp > Binv[i][m] / p_hat[i] )
			{
				temp = Binv[i][m] / p_hat[i];
				*r = i;
			}
		}
	}

	if( temp == INFTY )
	{
		*r = -1;
//		fprintf( stderr, "r = -1!!!!\n" );
//		exit( EXIT_FAILURE );
	}
}


void b_find_r( int m, int n, double **Binv, double *p_hat, int *r , int *basis )
{
	int i;
	double temp;

	temp = INFTY;
	for( i = 0; i < n; i++ )
	{
		if( p_hat[i] > EPS && Binv[i][m] > -EPS )
		{
			if(  temp == Binv[i][m] / p_hat[i] && basis[i] < basis[*r] )
			{
				temp = Binv[i][m] / p_hat[i];
				*r = i;
				continue;
			}
			if( temp > Binv[i][m] / p_hat[i] )
			{
				temp = Binv[i][m] / p_hat[i];
				*r = i;
			}
		}
	}

	if( temp == INFTY )
	{
		*r = -1;
//		fprintf( stderr, "r = -1!!!!\n" );
//		exit( EXIT_FAILURE );
	}
}


void lex_find_r( int m, int n, double **Binv, double *p_hat, int *r, int eps_key )
{
	int i, j;
	double temp;
	double q0;		
	double qp;		

	
	temp = INFTY;
	for( i = 0; i < m; i++ )
	{
		if( p_hat[i] > EPS && Binv[i][m] > -EPS )
		{
			if( eps_key > 0 )
			{
				if( fabs( temp - Binv[i][m] / p_hat[i] ) <= DBL_EPSILON * abs_max( temp, Binv[i][m] / p_hat[i] ) )
				{
					for( j = 0; j < m; j++ )
					{
						q0 = Binv[i][j] / p_hat[i];		
						qp = Binv[*r][j] / p_hat[*r];	
						if( q0 != qp )								 
							break;
					}
					if( q0 < qp )									
					{
						temp = Binv[i][m] / p_hat[i];
 						*r = i;
					}
					continue;
				}
			}
			else if( eps_key == 0 )
			{
				if( fabs( temp - Binv[i][m] / p_hat[i] ) < EPS )
				{
					for( j = 0; j < m; j++ )
					{
						q0 = Binv[i][j] / p_hat[i];		
						qp = Binv[*r][j] / p_hat[*r];	
						if( q0 != qp )								 
							break;
					}
					if( q0 < qp )									
					{
						temp = Binv[i][m] / p_hat[i];
 						*r = i;
					}
					continue;
				}
			}
			else
			{
				if( temp == Binv[i][m] / p_hat[i] )
				{
					for( j = 0; j < m; j++ )
					{
						q0 = Binv[i][j] / p_hat[i];		
						qp = Binv[*r][j] / p_hat[*r];	
						if( q0 != qp )								
							break;
					}
					if( q0 < qp )									
					{
						temp = Binv[i][m] / p_hat[i];
 						*r = i;
					}
					continue;
				}
			}
			if( temp > Binv[i][m] / p_hat[i] )
			{
				temp = Binv[i][m] / p_hat[i];
				*r = i;
			}
		}
	}

	if( temp == INFTY )
	{
		*r = -1;
//		fprintf( stderr, "r = -1!!!!\n" );
//		exit( EXIT_FAILURE );
	}
}



int degeneracy( int m, double **Binv )
{
	int i;

	for( i = 0; i < m; i++ )
	{
		if( Binv[i][m] == 0.0 ) return 1;		
//		if( ( Binv[i][m] >= -EPS ) && ( Binv[i][m] <= EPS ) ) return 1;
	}
	return 0;
}


int phase1( int m, int n, double **Binv, double **p, double *d, int *basis, int *basisset )
{
	int i, r, s;
	double d_s, *p_hat;
	p_hat = vec( m+1 );

	while( 1 )
	{	
		find_s( m+1, n, Binv, p, d, basisset, &s, &d_s );

		if( d_s >= -EPS )
		{
			free( p_hat );
			if( Binv[m+1][m+1] < -EPS ) return( 1 );
			else return( 0 );
		}
		for( i = 0; i < m + 1; i++ ) p_hat[i] = inner( Binv[i], p[s], m+1 );
		find_r( m+1, m, Binv, p_hat, &r );

		if( basis[r] < n ) basisset[basis[r]] = 0;
		basis[r] = s;
		basisset[s] = 1;

		pivoting( m+1, Binv, p_hat, d_s, r );
	}
}



/* Arrange phase1 to simply generate a valid nonbasis information*/
int u_phase1(int m, int n, double **Binv, double **p, double *d, int *basis, int *basisset)
{
	int i, r, s, bind, j, ukey;
	double d_s, *p_hat, *u;
	p_hat = vec(m + 1);

	bind = n - m; //number of binding equations (minus 1 for artificial var)
	while (1)
	{
		find_s(m + 1, n, Binv, p, d, basisset, &s, &d_s);
		if (d_s >= -EPS)
		{
			free(p_hat);
			if (Binv[m + 1][m + 1] < -EPS) return(-1);  //if phase1 fails.
			else{  //if phase1 passes, it gives a valid nonbasis info		
				u = vec(n);							// initialize u_vec solution
				//printvec(u, n);
				ukey = 0;							// ukey to represent a valid u.
				//printf("\n");
				//printmat(Binv, m+2,m+2);
				for (i = 0; i < m; i++){ u[basis[i]] = Binv[i][m + 1]; }
				//solution(m, nb, Binv, basis, u);	// retrieve solution for u
				for (j = 0; j < bind; j++){
					if (u[j] > EPS){ ukey += 1 << j; }
				}

				//printf("In the u_phase1, bind is %d,  basis is  ", bind);
				//i_printvec(basis, m+1);
				//printvec(u, n);
				//printf(", and ukey is %d \t", ukey);
				//pbinf(ukey);
				free(u);
				return(ukey);
			}
		}
		for (i = 0; i < m + 1; i++) p_hat[i] = inner(Binv[i], p[s], m + 1);
		find_r(m + 1, m, Binv, p_hat, &r);


		if (basis[r] < n) basisset[basis[r]] = 0;
		basis[r] = s;
		basisset[s] = 1;

		pivoting(m + 1, Binv, p_hat, d_s, r);
	}
}



int simple_phase1(int m, int n, double **A, double *b)
{
	int ukey; // register int i;
	double **p, **Binv, *c, *d;
	int *basis, *basisset;

	p = mat(n + m, m + 1);
	Binv = mat(m + 2, m + 2);
	c = vec(n + m);
	d = vec(n + m);
	basis = ivec(m);
	basisset = ivec(n + m);

	tableau(m, n, Binv, p, A, b, c, d, basis, basisset);
	ukey = u_phase1(m, n + m, Binv, p, d, basis, basisset);

	freemat(p, n + m);
	freemat(Binv, m + 2);
	free(c);
	free(d);
	free(basis);
	free(basisset);
	return(ukey);
}



int phase2( int m, int n, double **Binv, double **p, double *c, int *basis, int *basisset, int avoid )
{
	int i, r, s;

		int cycle = 0;

	
	int count_cycle = 0;
	double before, after;

	double c_s, *p_hat;


	p_hat = vec( m );

	while( 1 )
	{
		if( cycle == 1 && avoid > 0 )
		{
			b_find_s( m, n, Binv, p, c, basisset, &s, &c_s );
		}
		else
			find_s( m, n, Binv, p, c, basisset, &s, &c_s );


		if( c_s >= -EPS )
		{
			free( p_hat );
			return( 0 );
		}
    
	
		for( i = 0; i < m; i++ ) p_hat[i] = inner( Binv[i], p[s], m );

		if( maxvalue( p_hat, m ) <= EPS )
		{
			free( p_hat );
			return( 2 );
		}


		if( avoid > 0 )
		{
			if( cycle )
				b_find_r( m, m, Binv, p_hat, &r, basis );
			else
				find_r( m, m, Binv, p_hat, &r );
		}
		else
			lex_find_r( m, m, Binv, p_hat, &r, -1 );


		if( r == -1 )
		{
			free( p_hat );
			return( 2 );
		}


		before = Binv[m][m];

		if( basis[r] < n ) basisset[basis[r]] = 0;
		basis[r] = s;
		basisset[s] = 1;

		pivoting( m, Binv, p_hat, c_s, r );



		if( avoid > 0 )
		{
			after = Binv[m][m];

			if( after == before )
			{
				count_cycle++;		
			}
			else
			{
				count_cycle = 0;	
				cycle = 0;
			}

			
			if( count_cycle >= avoid )
			{
				count_cycle = 0;
				cycle = 1;
			}
		}

	}
}



int easy_phase1( int m, int n, double **A, double *b )
{
	register int i;
	double **p, **Binv, *c, *d;
	int *basis, *basisset;

	p = mat( n+m, m+1 );
	Binv = mat( m+2, m+2 );
	c = vec( n+m );
	d = vec( n+m );
	basis = ivec( m );
	basisset = ivec( n+m );

	tableau( m, n, Binv, p, A, b, c, d, basis, basisset );
	i = phase1( m, n+m, Binv, p, d, basis, basisset );
	freemat( p, n+m );
	freemat( Binv, m+2 );
	free( c );
	free( d );
	free( basis );
	free( basisset );
	return( i );
}

int simplex( int m, int n, double **Binv, double **p, double **A, double *b, double *c, int *basis, int *basisset, int avoid )
{
	int i,j;
	double **Binv2, *d;

	d = vec( n+m );
	Binv2 = mat( m+2, m+2 );
	tableau( m, n, Binv2, p, A, b, c, d, basis, basisset );

	if( phase1( m, m+n, Binv2, p, d, basis, basisset ) != 0 )
	{    
		freemat( Binv2, m+2 );
		free( d );
		return( 1 );
	}

	for( j = 0; j < n + m; j++ )
	{
		if( d[j] + inner( Binv2[m+1], p[j], m ) > EPS ) basisset[j] = 2;
	}
	for( i = 0; i < m + 1; i++)
		for( j = 0; j < m + 1; j++ ) Binv[i][j] = Binv2[i][j];
	for( i = 0; i < m + 1; i++ ) Binv[i][m] = Binv2[i][m+1];
	freemat( Binv2, m+2 );
	free( d );

	return( phase2( m, m+n, Binv, p, c, basis, basisset, avoid ) );
}

double easy_simplex( int m, int n, double **A, double *b, double *c, double *x, int avoid )
{
	int i, j;
	double **Binv, *d, temp, **p;
	int *basis, *basisset;

	d = vec( n+m );
	Binv = mat( m+2, m+2 );
	p = mat( n+m, m+1 );
	basis = ivec( m );
	basisset = ivec( n+m );

	tableau( m, n, Binv, p, A, b, c, d, basis, basisset );
	if( phase1( m, m+n, Binv, p, d, basis, basisset ) == 1 )
	{
		printf( "!!!!\n" );
		exit( 1 );
	}
	for( j = 0; j < m + n; j++ )
	{
		if( d[j] + inner( Binv[m+1], p[j], m ) > EPS ) basisset[j] = 2;
	}
	for( i = 0; i < m + 1; i++ ) Binv[i][m] = Binv[i][m+1];
	phase2( m, m+n, Binv, p, c, basis, basisset, avoid );
	for( i = 0; i < n; i++ ) x[i] = 0;
	for( i = 0; i < m; i++ )
		if( basis[i] < n ) x[basis[i]] = Binv[i][m];
	temp = -Binv[m][m];
	freemat( Binv, m+2 );
	freemat( p, n+m );
	free( d );
	free( basis );
	free( basisset );
	return( temp );
}

double r_easy_simplex( int m, int n, double **A, double *b, double *c, int avoid )
{
	int i, j;
	double **Binv, *d, temp, **p;
	int *basis, *basisset;

	d = vec( n+m );
	Binv = mat( m+2, m+2 );
	p = mat( n+m, m+1 );
	basis = ivec( m );
	basisset = ivec( n+m );

	tableau( m, n, Binv, p, A, b, c, d, basis, basisset );
	if( phase1( m, m+n, Binv, p, d, basis, basisset ) == 1 )
	{
		printf( "!!!!\n" );
		exit( 1 );
	}
	for( j = 0; j < m + n; j++ )
	{
		if( d[j] + inner( Binv[m+1], p[j], m ) > EPS ) basisset[j] = 2;
	}
	for( i = 0; i < m + 1; i++ ) Binv[i][m] = Binv[i][m+1];

	if( phase2( m, m+n, Binv, p, c, basis, basisset, avoid ) == 2 )
		temp = LARGE;
	else
		temp = Binv[m][m];

	freemat( Binv, m+2 );
	freemat( p, n+m );
	free( d );
	free( basis );
	free( basisset );
	return( temp );
}

