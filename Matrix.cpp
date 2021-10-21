#include "matrix.h"
#include <cstring>
#include <cmath>
#include <vector>
#include <cstdarg>
#include "exceptions.h"
#include "algorithm"

using std::vector;

custom::matrix::matrix( unsigned int m /*= 0*/, unsigned int n /*= 0*/ )
{
	this->data = NULL;
	this->eps = 1e-10;
	this->rank = 0;

	if ((0 == m)||(0 == n)) 
	{
		this->lines = this->columns = 0;
		return;
	}

	this->lines = m;
    this->columns = n;

	this->data = new double[this->lines*this->columns];
	memset(this->data,0,this->lines*this->columns*sizeof(double)); /* fill with zero */
}

custom::matrix::matrix( const custom::matrix& mat )
{
	this->lines = mat.lines;
	this->columns = mat.columns;
	this->eps = mat.eps;
	this->rank = mat.rank;

	if ((0 == this->lines)||(0 == this->columns)) 
	{
		this->lines = this->columns = 0;
		this->data = NULL;
		return;
	}
	this->data = new double[this->lines*this->columns];
	memcpy(this->data,mat.data,this->lines*this->columns*sizeof(double));
}

custom::matrix::~matrix()
{
	clear();
}

void custom::matrix::resize( unsigned int m, unsigned int n )
{
	size_t cpr,cpc,oldr,oldc;
	oldr = this->lines;
    oldc = this->columns;
	cpr = (m < oldr) ? m : oldr;/* 需复制的元素行数 */
	cpc = (n < oldc) ? n : oldc;/* 需复制的列数 */

	double *p = new double[m*n];
	/* copy the original data */
	for (size_t i = 0; i < cpr; i++)
	{
		for (size_t j = 0; j < cpc; j++)
		{
			p[i*n+j] = data[i*oldc+j];
		}
	}
	/* delete the original memory */
	this->clear();
	/* set data */
	this->data = p;
	this->lines = m;
	this->columns = n;
}

void custom::matrix::clear()
{
	if (data != NULL) delete[] data;
	data = NULL;
	this->lines = 0;
	this->columns = 0;
}

double& custom::matrix::at( unsigned int i, unsigned int j )
{
	if (i >= this->lines) {throw out_of_range_exception();}
	if (j >= this->columns) {throw out_of_range_exception();}
	return data[i*this->columns+j];
}

double custom::matrix::at( unsigned i, unsigned j ) const
{
	if (i >= this->lines) {throw out_of_range_exception();}
	if (j >= this->columns) {throw out_of_range_exception();}
	return data[i*this->columns+j];
}

/*
double* custom::matrix::operator[]( unsigned int i )
{
	return data[i];
}*/

custom::matrix custom::matrix::operator*( const custom::matrix& A ) const
{
    custom::matrix M;
	if (this->columns != A.lines)
	{
		M.clear();
        throw generic_exception("matrix dimension not match!");
	}
	else
	{
		M.resize(this->lines, A.columns);
		for (unsigned int i = 0; i < M.get_lines(); i++)
		{
			for (unsigned int j = 0; j < M.get_columns(); j++)
			{
				M[i][j] = 0;
				for (unsigned int k = 0; k < this->columns; k++)
				{
					M[i][j] += data[i*this->columns+k] * A.at(k,j);
				}
			}
		}
	}
	return M;
}

custom::matrix custom::matrix::operator*( double a )
{
    custom::matrix r(this->lines,this->columns);
	for (unsigned int i = 0; i < this->lines; i++)
	{
		for (unsigned int j = 0; j < this->columns; j++)
		{
			r[i][j] = this->data[i*this->columns+j]*a;
		}
	}
	return r;
}

custom::matrix& custom::matrix::operator=(const custom::matrix& mat)
{
	clear();
	this->lines = mat.lines;
	this->columns = mat.columns;
	this->eps = mat.eps;
	this->rank = mat.rank;

	this->data = new double[this->lines*this->columns];
	memcpy(this->data, mat.data, this->lines*this->columns*sizeof(double));

	return *this;
}

custom::matrix& custom::matrix::operator=(double* data)
{
	memcpy(this->data, data, this->lines*this->columns*sizeof(double));
	return *this;
}

custom::matrix custom::matrix::inverse()
{
    if(this->lines != this->columns) return custom::matrix(0,0);

	double *mat = new double[lines*columns];
	unsigned int *maprow = new unsigned int[lines];
    // prepare the matrix and the mapping
	for (unsigned int i = 0; i < lines; i++)
	{
		maprow[i] = i;
		for (unsigned int j = 0; j < columns; j++)
		{
			mat[i*columns+j] = data[i*columns+j];
		}
	}
	/* begin the actual process, using the Gauss-Jordan method */
	double pivot;
	unsigned int pivot_index;
	for (unsigned int i = 0; i < lines; i++)
	{
		/* find the pivot of the column */
		pivot = mat[maprow[i]*columns+i];
		pivot_index = i;
		for (unsigned int j = i; j < lines; j++)
		{
			if(fabs(mat[maprow[j]*columns+i]) > fabs(pivot)) {pivot = mat[maprow[j]*columns+i];pivot_index = j;} 
		}
		/* swap the two lines, i and pivot_index */
		unsigned int tmp = maprow[i];
		maprow[i] = maprow[pivot_index];
		maprow[pivot_index] = tmp;
		/* calculate the ith line */
		for (unsigned int j = 0; j < columns; j++)
		{
			if(j != i) mat[maprow[i]*columns+j] /= mat[maprow[i]*columns+i];
		}
		mat[maprow[i]*columns+i] = 1.0/mat[maprow[i]*columns+i];
		/* calculate other lines */
		for (unsigned int k = 0; k < lines; k++)
		{
			if(k == i) continue;
			for (unsigned int l = 0; l < columns; l++)
			{
				if(l == i) continue;
				mat[maprow[k]*columns+l] -= mat[maprow[i]*columns+l]*mat[maprow[k]*columns+i];
			}
			mat[maprow[k]*columns+i] = -mat[maprow[k]*columns+i]*mat[maprow[i]*columns+i];
		}
	}
	/* rearrange the columns */
    custom::matrix ivs(lines,columns);
	for (unsigned int i = 0; i < lines; i++)
	{
		for (unsigned int j = 0; j < columns; j++)
		{
			ivs[i][j] = mat[maprow[i]*columns+maprow[j]];
		}
	}

	delete[] maprow;
	delete[] mat;
	return ivs;
}

custom::matrix custom::matrix::inverse4X4()
{
    custom::matrix inverseMatrix(3,3);
    for(size_t i = 0; i < 3; ++i)
    {
        inverseMatrix.setValue(i, 0, getValue(i, 0));
        inverseMatrix.setValue(i, 1, getValue(i, 1));
        inverseMatrix.setValue(i, 2, getValue(i, 2));
    }
    custom::matrix inverseMatrix1 = inverseMatrix.transpose();
    custom::matrix reInverseMatrix(4,4);
    for(size_t i = 0; i < 3; ++i)
    {
        reInverseMatrix.setValue(i, 0, inverseMatrix1.getValue(i, 0));
        reInverseMatrix.setValue(i, 1, inverseMatrix1.getValue(i, 1));
        reInverseMatrix.setValue(i, 2, inverseMatrix1.getValue(i, 2));
        reInverseMatrix.setValue(i, 3, 0);
    }
    // 下面求 原点在三个坐标轴的投影
    float dx = getValue(3,0)*getValue(0,0) + getValue(3,1)*getValue(0,1) + getValue(3,2)*getValue(0,2);
    float dy = getValue(3,0)*getValue(1,0) + getValue(3,1)*getValue(1,1) + getValue(3,2)*getValue(1,2);
    float dz = getValue(3,0)*getValue(2,0) + getValue(3,1)*getValue(2,1) + getValue(3,2)*getValue(2,2);
    reInverseMatrix.setValue(3, 0, -dx);
    reInverseMatrix.setValue(3, 1, -dy);
    reInverseMatrix.setValue(3, 2, -dz);
    reInverseMatrix.setValue(3, 3, 1);
    return reInverseMatrix;
}

custom::matrix custom::matrix::transpose()
{
    custom::matrix T(this->columns,this->lines);
	for (unsigned int i = 0; i < columns; i++)
	{
		for (unsigned int j = 0; j < lines; j++)
		{
			T[i][j] = data[j*columns+i];
		}
	}
	return T;
}

custom::matrix custom::matrix::Hermite_canonical_form()
{
	double *mat = new double[lines*columns];
	unsigned int *maprow = new unsigned int[lines];
    /* prepare the matrix for fast access */
	for (unsigned int i = 0; i < lines; i++)
	{
		maprow[i] = i;
		for (unsigned int j = 0; j < columns; j++)
		{
			mat[i*columns+j] = this->at(i,j);
		}
	}
	/* 初等变换化为标准型 */
	double pivot;
	unsigned int pivot_index;
	for (unsigned int i = 0; i < lines; i++)
	{
		/* find first non-zero column */
		int nzcol = -1;// indicate the non-zero column
		for (unsigned int j = i; j < columns; j++)
		{
			for (unsigned int k = i; k < lines; k++)
			{
				if(fabs(mat[maprow[k]*columns+j]) > this->eps) {nzcol = j;break;} /// @todo 判断浮点数非零,小心
			}
			if(nzcol >= 0) break;
		}
		if (nzcol < 0) // can't find the non-zero column ,finish the process
		{
			this->rank = i;
			break;
		}
		/* find the pivot of the line */
		pivot = mat[maprow[i]*columns+nzcol];
		pivot_index = i;
		for (unsigned int j = i; j < lines; j++)
		{
			if(fabs(mat[maprow[j]*columns+nzcol]) > fabs(pivot)) {pivot = mat[maprow[j]*columns+nzcol];pivot_index = j;} 
		}
		/* swap the two lines, i and pivot_index */
		unsigned int tmp = maprow[i];
		maprow[i] = maprow[pivot_index];
		maprow[pivot_index] = tmp;
		/* calculate the ith line */
		for (unsigned int j = nzcol+1; j < columns; j++)
		{
			mat[maprow[i]*columns+j] /= mat[maprow[i]*columns+nzcol];
		}
		mat[maprow[i]*columns+nzcol] = 1.0;
		/* calculate other lines */
		for (unsigned int k = 0; k < lines; k++)
		{
			if(k == i) continue;
			for (unsigned int j = nzcol+1; j < columns; j++)
			{
				mat[maprow[k]*columns+j] -= mat[maprow[i]*columns+j]*mat[maprow[k]*columns+nzcol];
			}
			mat[maprow[k]*columns+nzcol] = 0;
		}
	}
    /* copy the custom::matrix */
    custom::matrix Hermite(lines,columns);
	for (unsigned int i = 0; i < lines; i++)
	{
		for (unsigned int j = 0; j < columns; j++)
		{
			Hermite[i][j] = mat[maprow[i]*columns+j];
		}
	}
	Hermite.rank = this->rank;

	delete[] maprow;
	delete[] mat;
	return Hermite;
}

custom::matrix custom::matrix::pseudo_inverse()
{
    custom::matrix L,R;
    custom::matrix::FullRankDecomposition(*this,L,R);
    custom::matrix LT,RT;
	/* L,R的转置 */
	LT = L.transpose();
	RT = R.transpose();
	/* 相乘 */
	L = LT*L;
	R = R*RT;
	/* 求逆 */
	L = L.inverse();
	R = R.inverse();
	/* B的广义逆 */
    custom::matrix Bplus;
	Bplus = ((RT*R)*L)*LT;
	return Bplus;
}

bool custom::matrix::null() {return (lines == 0)||(columns == 0);}

void custom::matrix::set_precision(double pc) {this->eps = pc;}

void custom::matrix::FullRankDecomposition( custom::matrix& B, custom::matrix& L, custom::matrix& R )
{
	//throw std::exception("The method or operation is not implemented.");;
    custom::matrix HB;
	HB = B.Hermite_canonical_form();
	vector<unsigned int> ones;
	/* find position of ones */
	for (unsigned int i = 0; i < HB.get_lines(); i++)
	{
		unsigned int j = 0;
		for (j = 0; j < HB.get_columns(); j++)
		{
			if(HB[i][j] == 1.0) 
			{
				ones.push_back(j);
				break;
			}
		}
		if(j == HB.get_columns()) break;
	}
	/* construct the L and R where B=LR */
    /* matrix L */
	L.resize(B.get_lines(),ones.size());
	for (unsigned int i = 0; i < L.get_lines(); i++)
	{
		for (unsigned int j = 0; j < L.get_columns(); j++)
		{
			L[i][j] = B[i][ones[j]];
		}
	}
    /* matrix R */
	R.resize(ones.size(),B.get_columns());
	for (unsigned int i = 0; i < R.get_lines(); i++)
	{
		memcpy(R[i],HB[i],R.get_columns()*sizeof(double));// copy lines
	}
}

custom::matrix custom::matrix::operator+( custom::matrix& mat )
{
    if((this->columns != mat.columns)||(this->lines != mat.lines)) throw generic_exception("operator+: matrix not the same size");
    custom::matrix r(this->lines,this->columns);
	for (unsigned int i = 0; i < r.get_lines(); i++)
	{
		for (unsigned int j = 0; j < r.get_columns(); j++)
		{
			r[i][j] = this->at(i,j) + mat[i][j];
		}
	}
	return r;
}

custom::matrix& custom::matrix::operator+=( custom::matrix& mat )
{
    if((this->columns != mat.columns)||(this->lines != mat.lines)) throw generic_exception("operator+=: matrix not the same size");
	for (unsigned int i = 0; i < this->lines; i++)
	{
		for (unsigned int j = 0; j < this->columns; j++)
		{
			this->at(i,j) += mat[i][j];
		}
	}
	return *this;
}

custom::matrix custom::matrix::operator-( custom::matrix& mat )
{
    if((this->columns != mat.columns)||(this->lines != mat.lines)) throw generic_exception("operator-: matrix not the same size");
    custom::matrix r(this->lines,this->columns);
	for (unsigned int i = 0; i < r.get_lines(); i++)
	{
		for (unsigned int j = 0; j < r.get_columns(); j++)
		{
			r[i][j] = this->at(i,j) - mat[i][j];
		}
	}
	return r;
}

custom::matrix custom::matrix::operator-()
{
    custom::matrix r(this->lines,this->columns);
	for (unsigned int i = 0; i < r.lines; i++)
	{
		for (unsigned int j = 0; j < r.columns; j++)
		{
			r[i][j] = -this->at(i,j);
		}
	}
	return r;
}

custom::matrix& custom::matrix::operator-=( custom::matrix& mat )
{
    if((this->columns != mat.columns)||(this->lines != mat.lines)) throw generic_exception("operator-=: matrix not the same size");
	for (unsigned int i = 0; i < this->lines; i++)
	{
		for (unsigned int j = 0; j < this->columns; j++)
		{
			this->at(i,j) -= mat[i][j];
		}
	}
	return *this;
}

void custom::matrix::init( double d, ... )
{
	std::va_list args;
	double ls = d;
	va_start(args, d);
		for (unsigned int i = 0; i < this->lines; i++)
		{
			for (unsigned int j = 0; j < this->columns; j++)
			{
				this->at(i,j) = ls;
				ls = va_arg(args,double);
			}
		}
	va_end(args);
}

unsigned int custom::matrix::get_columns() const {return this->columns;}

unsigned int custom::matrix::get_lines() const {return this->lines;}

const double *custom::matrix::get_data() const {return this->data;}

double& custom::matrix::operator()( unsigned int i, unsigned int j )
{
    if (i >= this->lines) {throw out_of_range_exception();}
	if (j >= this->columns) {throw out_of_range_exception();}
	return data[i*this->columns+j];
}

double *custom::matrix::operator[](unsigned int i) {return &data[i*this->columns];}

void custom::matrix::load_identity()
{
    int sz = std::min(this->lines,this->columns);
	memset(data,0,this->lines*this->columns*sizeof(double));
	for (int i = 0; i < sz; i++)
	{
		at(i,i) = 1;
	}
}

std::ostream& operator<<(std::ostream& out, custom::matrix& mat)
{
	for (unsigned int i = 0; i < mat.get_lines(); i++)
	{
		for (unsigned int j = 0; j < mat.get_columns(); j++)
		{
			out<<mat[i][j]<<'\t';
		}
		out<<std::endl;
	}
	return out;
}

void  custom::matrix::setrowVal(int row, float vrow[])
{
    for(int j=0;j<(int)this->columns;j++)
    {
        data[row*this->columns+j] = vrow[j];
    }
}

void  custom::matrix::setValue(int row, int col, float val)
{
    data[row*columns+col] = val;
}

float custom::matrix::getValue(int row, int col)
{
    return data[row*columns+col];
}
