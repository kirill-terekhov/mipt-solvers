#ifndef _CSRMATRIX_H
#define _CSRMATRIX_H

#include "datatype.h"
#include "vector.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <limits>


template<typename KeyType>
class CSRMatrixType
{
	std::vector<idx_t> ia, ja;
	std::vector<KeyType> a;
public:
	CSRMatrixType() {ia.push_back(0);}
	CSRMatrixType(const std::vector<KeyType> & diag)
	{
		ia.resize(diag.size()+1);
		ja.resize(diag.size());
		a.resize(diag.size());
		for(idx_t k = 0; k < (idx_t)diag.size(); ++k)
		{
			a[k] = diag[k];
			ja[k] = k;
			ia[k] = k;
		}
		ia.back() = (idx_t)diag.size();
	}
	CSRMatrixType(const std::vector<idx_t> & ia, const std::vector<idx_t> & ja, const std::vector<KeyType> & a)
	: ia(ia), ja(ja), a(a) {}
	CSRMatrixType(const CSRMatrixType & b) : ia(b.ia), ja(b.ja), a(b.a) {}
	CSRMatrixType & operator=(CSRMatrixType const & b) { ia = b.ia; ja = b.ja; a = b.a; return *this;}
	void Clear()
	{
		ia.clear();
		ja.clear();
		a.clear();
		ia.push_back(0);
	}
	void LoadBinary(std::string file)
	{
		std::ifstream input(file.c_str(), std::ios::binary);
		size_t header[3];
		input.read(reinterpret_cast<char *>(header),sizeof(size_t)*3);
		ia.resize(header[0]);
		ja.resize(header[1]);
		a.resize(header[2]);
		input.read(reinterpret_cast<char *>(&ia[0]),sizeof(idx_t)*ia.size());
		input.read(reinterpret_cast<char *>(&ja[0]),sizeof(idx_t)*ja.size());
		input.read(reinterpret_cast<char *>(&a[0]),sizeof(KeyType)*a.size());
		input.close();
	}
	void SaveBinary(std::string file) const
	{
		std::ofstream output(file.c_str(), std::ios::binary);
		size_t header[3] = {ia.size(),ja.size(),a.size()};
		output.write(reinterpret_cast<const char *>(header),sizeof(size_t)*3);
		output.write(reinterpret_cast<const char *>(&ia[0]),sizeof(idx_t)*ia.size());
		output.write(reinterpret_cast<const char *>(&ja[0]),sizeof(idx_t)*ja.size());
		output.write(reinterpret_cast<const char *>(&a[0]),sizeof(KeyType)*a.size());
		output.close();
	}
	void Load(std::string file)
	{
		idx_t N,M,L;
		idx_t row, col, last_row;
		KeyType val;
		std::ifstream input(file.c_str());
		if( input.fail() )
		{
			std::cout << "Cannot open " << file << std::endl;
			exit(-1);
		}
		ia.clear();
		ja.clear();
		a.clear();
		while( !input.eof() )
		{
			int c = input.peek();
			if( isspace(c) ) 
			{
				c = input.get();
				continue;
			}
			if( c == '%' ) 
			{
				input.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
				continue;
			}
			input >> N >> M >> L;
			//assert(N == M);
			ia.resize(N+1);
			ja.reserve(L);
			a.reserve(L);
			break;
		}
		std::cout << "Sizes: " << N << " " << M << " nonzeros " << L << std::endl;
		if( ia.empty() ) 
			throw "MatrixMarket header not found";
		last_row = 0;
		ia[0] = 0;
		while( !input.eof() )
		{
			int c = input.peek();
			if( isspace(c) ) 
			{
				c = input.get();
				continue;
			}
			if( c == '%' )
			{
				input.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
				continue;
			}
			if (input.eof())
				break;
			input >> row >> col >> val;
			//if (input.eof() || input.fail())
			//	break;
			if( !(row == last_row || row == last_row+1) )
				throw "Matrix rows are unordered";
			if( row != last_row )
			{
				ia[row] = ia[last_row];
				last_row = row;
			}
			ia[row]++;
			ja.push_back(col-1);
			a.push_back(val);
		}
		input.close();
	}
	void Save(std::string file) const
	{
		std::ofstream output(file.c_str());
		output << "%%MatrixMarket matrix coordinate real general" << std::endl;
		output << Size() << " " << Size() << " " << ja.size() << std::endl;
		output << std::scientific;
		output.precision(16);
		for(idx_t i = 0; i < Size(); ++i)
		{
			for(idx_t j = ia[i]; j < ia[i+1]; ++j)
				output << i+1 << " " << ja[j]+1 << " " << a[j] << std::endl;
		}
		output.close();
	}
	// y = alpha*A*x + beta*y
	void Multiply(const KeyType & alpha, const std::vector<KeyType> & x, const KeyType & beta, std::vector<KeyType> & y) const
	{
		if(y.size() != Size()) throw "Wrong argument size";
		//TODO: check columns in A equal to size in x
		for(idx_t i = 0; i < static_cast<idx_t>(Size()); ++i)
		{
			y[i] *= beta;
			for(idx_t j = ia[i]; j < ia[i+1]; ++j)
				y[i] += alpha*a[j]*x[ja[j]];
		}
	}
	// y = alpha*A.Transpose()*x + beta*y
	void MultiplyTranspose(const KeyType & alpha, const std::vector<KeyType> & x, const KeyType & beta, std::vector<KeyType> & y) const
	{
		if(x.size() != Size()) throw "Wrong argument size";
		for(idx_t i = 0; i < static_cast<idx_t>(y.size()); ++i)
			y[i] *= beta;
		for(idx_t i = 0; i < static_cast<idx_t>(Size()); ++i)
		{
			for(idx_t j = ia[i]; j < ia[i+1]; ++j)
				y[ja[j]] += alpha*a[j]*x[i];
		}
	}
	idx_t Columns(bool square = false) const
	{
		idx_t Cols = 0;
		if( square )
			Cols = Size();
		else
		{
			idx_t Colsloc = 0;
			for(idx_t i = 0; i < static_cast<idx_t>(ja.size()); ++i) 
				Colsloc = std::max(Colsloc,ja[i]+1);
			Cols = std::max(Cols,Colsloc);
		}
		return Cols;
	}
	CSRMatrixType Transpose(bool square = false) const
	{
		idx_t Cols = Columns(square);
		std::vector< idx_t > iat(Cols+1), ist(Cols,0);
		std::vector< idx_t > jat(Nonzeros()), pos(Nonzeros());
		std::vector<KeyType> at(a.size());
		//find out row size
		for(idx_t i = 0; i < Nonzeros(); ++i)
			pos[i] = ist[ja[i]]++;
		//setup ia array and set writing position
		iat[0] = 0;
		for(idx_t i = 0; i < Cols; ++i)
			iat[i+1] = iat[i] + ist[i];
		for(idx_t i = 0; i < static_cast<idx_t>(Size()); ++i)
		{
			for(idx_t j = ia[i]; j < ia[i+1]; ++j)
				jat[iat[ja[j]] + pos[j]] = i;
		}
		for(idx_t i = 0; i < static_cast<idx_t>(Nonzeros()); ++i)
			at[iat[ja[i]] + pos[i]] = a[i];
		return CSRMatrixType(iat,jat,at);
	}
	void Print(int w = 5) const
	{
		idx_t Cols = Columns();
		for(idx_t i = 0; i < Size(); ++i)
		{
			std::cout << "|";
			idx_t k = 0;
			for(idx_t j = ia[i]; j < ia[i+1]; ++j)
			{
				while(k++ < ja[j]) std::cout << std::setw(w) << " ";
				std::cout << std::setw(w) << a[j];
			}
			while(k++ < Cols) std::cout << std::setw(w) << " ";
			std::cout << "|" << std::endl;
		}
	}
	void PrintContents() const
	{
		std::cout << "ia: ";
		for(idx_t i = 0; i < (idx_t)ia.size(); ++i) std::cout << ia[i] << " ";
		std::cout << std::endl;
		std::cout << "ja: ";
		for(idx_t i = 0; i < (idx_t)ja.size(); ++i) std::cout << ja[i] << " ";
		std::cout << std::endl;
		std::cout << "a: ";
		for(idx_t i = 0; i < (idx_t)a.size(); ++i) std::cout << a[i] <<  " ";
		std::cout << std::endl;
	}
	void Diagonal(std::vector<KeyType> & d) const
	{
		d.resize(Size());
		for(idx_t i = 0; i < static_cast<idx_t>(Size()); ++i)
		{
			d[i] = 0;
			for(idx_t j = 0; j < RowSize(i); ++j) if( Col(i,j) == i ) 
			{
				d[i] = Val(i,j); 
				break;
			}
		}
	}
	bool DiagonalPosition(std::vector<idx_t> & d) const
	{
		d.resize(Size());
		idx_t err = 0;
		for(idx_t i = 0; i < Size(); ++i)
		{
			d[i] = RowSize(i);
			for(idx_t j = 0; j < RowSize(i); ++j) if( Col(i,j) == i ) 
			{
				d[i] = j; 
				break;
			}
			if( d[i] == RowSize(i) )
				err++;
		}
		return err == 0;
	}
	void SortRows()
	{
		for(idx_t k = 0; k < Size(); ++k)
		{
			idx_t n = ia[k+1] - ia[k], s = ia[k];
			if( n )
			{
				for (idx_t i = 0; i < n-1; i++)
				{
					for (idx_t j = 0; j < n-i-1; j++)
					{
						if (ja[s+j] > ja[s+j+1])
						{
							std::swap(ja[s+j],ja[s+j+1]);
							std::swap( a[s+j], a[s+j+1]);
						}
					}
				}
			}
		}
	}
	/// Check if the matrix is diagonal for a given tolernace.
	bool Diagonal(double eps = 1.0e-7) const
	{
		bool ret = true;
		std::vector<KeyType> diag;
		Diagonal(diag);
		for(idx_t i = 0; i < Size(); ++i)
			diag[i] = std::fabs(diag[i]);
		for(idx_t i = 0; i < Size(); ++i)
		{
			for(idx_t j = 0; j < RowSize(i); ++j)
				if( Col(i,j) != i )
					ret &= fabs(Val(i,j)) < eps*std::max(diag[i],diag[Col(i,j)]);
		}
		return ret;
	}
	/// Check if the matrix is triangular for a given tolernace.
	/// Diagnal matrix is also considered triangular.
	bool Triangular(double eps = 1.0e-7) const
	{
		KeyType Lnorm = 0, Unorm = 0, Dnorm = 0;
		for(idx_t i = 0; i < Size(); ++i)
		{
			for(idx_t j = 0; j < RowSize(i); ++j)
				if( Col(i,j) < i )
					Lnorm += std::fabs(Val(i,j));
				else if( Col(i,j) > i )
					Unorm += std::fabs(Val(i,j));
				else
					Dnorm += std::fabs(Val(i,j));
		}
		if( Unorm > eps * Dnorm && Lnorm > eps * Dnorm ) // both upper and lower triangular are nonzero
			return false;
		return true;
	}
	/// Check if the matrix is symmetric for a given tolernace.
	bool Symmetric(double eps = 1.0e-7) const
	{
		KeyType normp = 0, normm = 0;
		std::vector<bool> mrk(Nonzeros(),false);
		for (idx_t i = 0; i < Size(); i++)
		{
			for (idx_t j = ia[i]; j < ia[i+1]; j++)
			{
				const idx_t i1 = ja[j];
				if( mrk[j] )
				{
					mrk[j] = !mrk[j];
					continue;
				}
				if (i1 < i)
				{
					normm += a[j] * a[j];
					normp += a[j] * a[j];
				}
				else if(i1 == i)
					normp += a[j] * a[j];
				else
				{
					for (idx_t j1 = ia[i1]; j1 < ia[i1+1]; j1++) if (ja[j1] == i)
					{
						mrk[j1] = !mrk[j1];
						normm += (a[j] - a[j1])*(a[j] - a[j1]);
						normp += (a[j] + a[j1])*(a[j] + a[j1]);
					}
				}
			}
		}
		std::cout << "sym: " << sqrt(normm/normp) << " m " << normm << " p " << normp << std::endl;
		if( normm < eps*eps*normp )
			return true;
		return false;
	}
	bool Compare(CSRMatrixType const & b, double eps = 1.0e-7) const
	{
		idx_t err = 0;
		KeyType diftol = 0;
		if( Size() == b.Size() )
		{
			{
				KeyType diftolloc = 0;
				for(idx_t i = 0; i < Size(); ++i)
				{
					KeyType aval = 0, bval = 0;
					for(idx_t j = ia[i]; j < ia[i+1]; ++j)
						aval += std::abs(a[j]);
					for(idx_t j = b.ia[i]; j < b.ia[i+1]; ++j)
						bval += std::abs(b.a[j]);
					KeyType div = std::max(std::numeric_limits<KeyType>::min(),std::abs(aval)+std::abs(bval));
					KeyType dif = std::abs(aval - bval);
					diftolloc = std::max(diftolloc,dif/(KeyType)div);
					if( dif > eps*div ) 
					{
						std::cout << "at row "<<i<< " row norms are " << aval << " and " << bval << " relative difference " << dif/div << std::endl;
						err++;
					}
					//~ else
					{
						for(idx_t j = ia[i]; j < ia[i+1]; ++j)
						{
							aval = a[j];
							bval = 0;
							for(idx_t j1 = b.ia[i]; j1 < b.ia[i+1]; ++j1)	if( ja[j] == b.ja[j1] )
							{
								bval = b.a[j1];
								break;	
							}
							KeyType div = std::max(std::numeric_limits<KeyType>::min(),std::abs(aval)+std::abs(bval));
							KeyType dif = std::abs(aval - bval);
							diftolloc = std::max(diftolloc,dif/div);
							if( dif > eps*div ) 
							{
								std::cout << "at A("<<i<<","<<ja[j]<<"), values are " << aval << " and " << bval << " relative difference " << dif/div << std::endl;
								err++;
							}
						}
					}
					{
						for(idx_t j = b.ia[i]; j < b.ia[i+1]; ++j)
						{
							aval = 0;
							bval = b.a[j];
							for(idx_t j1 = ia[i]; j1 < ia[i+1]; ++j1)	if( b.ja[j] == ja[j1] )
							{
								aval = a[j1];
								break;	
							}
							KeyType div = std::max(std::numeric_limits<KeyType>::min(),std::abs(aval)+std::abs(bval));
							KeyType dif = std::abs(aval - bval);
							diftolloc = std::max(diftolloc,dif/div);
							if( dif > eps*div ) 
							{
								std::cout << "at B("<<i<<","<<b.ja[j]<<"), values are " << aval << " and " << bval << " relative difference " << dif/div << std::endl;
								err++;
							}
						}
					}
				}
				diftol = std::max(diftol,diftolloc);
			}
			if( err == 0 ) 
				std::cout << "Matrices are same" << std::endl;
			else
				std::cout << "Relative difference tolerance " << diftol << std::endl;
		}
		else 
		{
			std::cout << "Matrix sizes " << Size() << " and " << b.Size() << std::endl;
			err++;
		}
		return err == 0;
	}
	inline size_t Nonzeros() const {return a.size();}
	inline idx_t Size() const {return (idx_t)(ia.size()-1);}
	inline idx_t RowSize(idx_t i) const { idx_t k = i + 1;  return ia[k] - ia[i]; }
	inline idx_t Col(idx_t i, idx_t k) const { idx_t q = ia[i] + k;  return ja[q]; }
	inline idx_t Col(const std::pair<idx_t,idx_t> & p) const {return Col(p.first,p.second);}
	inline idx_t & Col(idx_t i, idx_t k) {return ja[ia[i]+k];}
	inline idx_t & Col(const std::pair<idx_t,idx_t> & p) {return Col(p.first,p.second);}
	inline const KeyType& Val(idx_t i, idx_t k) const { idx_t q = ia[i] + k;  return a[q]; }
	inline const KeyType & Val(const std::pair<idx_t,idx_t> & p) const {return Val(p.first,p.second);}
	inline KeyType & Val(idx_t i, idx_t k) {return a[ia[i]+k];}
	inline KeyType & Val(const std::pair<idx_t,idx_t> & p) {return Val(p.first,p.second);}
	size_t Bytes() const {return get_bytes(ia) + get_bytes(ja) + get_bytes(a);}
	void Insert(const std::vector<idx_t> & ja_part, const std::vector<KeyType> & a_part) 
	{
		ja.insert(ja.end(),ja_part.begin(),ja_part.end());
		a.insert(a.end(),a_part.begin(),a_part.end());
	}
	void PushBack(idx_t col, const KeyType & val) 
	{ 
		ja.push_back(col);
		a.push_back(val);
	}
	void FinalizeRow() { ia.push_back((idx_t)ja.size()); }
	
	void ReserveSize(size_t size) {ia.reserve(size+1);}
	void ReserveNonzeros(size_t nonzeros) {ja.reserve(nonzeros); a.reserve(nonzeros);}
	void Shrink() { ia.shrink_to_fit(); ja.shrink_to_fit(); a.shrink_to_fit(); }
	
	inline std::vector<idx_t> & get_ia() {return ia;}
	inline std::vector<idx_t> & get_ja() {return ja;}
	inline std::vector<KeyType> & get_a() {return a;}
	inline const std::vector<idx_t> & get_ia() const {return ia;}
	inline const std::vector<idx_t> & get_ja() const {return ja;}
	inline const std::vector<KeyType> & get_a() const {return a;}
};

typedef CSRMatrixType<double> CSRMatrix;

double Resid(const CSRMatrix & A, const std::vector<double> & b, const std::vector<double> & x)
{
	double r = 0;
	for(idx_t i = 0; i < A.Size(); ++i)
	{
		double c = b[i];
		for(idx_t j = 0; j < A.RowSize(i); ++j)
			c -= A.Val(i,j)*x[A.Col(i,j)];
		r += c*c;
	}
	return std::sqrt(r);
}


#endif //_CSRMATRIX_H
