#ifndef _GRAPH_H
#define _GRAPH_H

#include "datatype.h"
#include "row_accumulator_graph.h"
#include <fstream>
#include <iostream>
#include <iomanip>

class CSRGraph
{
	std::vector<idx_t> ia, ja;
	std::string tmpfile;
public:
	CSRGraph() :ia(), ja(), tmpfile("") {ia.push_back(0);}
	CSRGraph(const std::vector<idx_t> & ia, const std::vector<idx_t> & ja)
	: ia(ia), ja(ja), tmpfile("") {}
	CSRGraph(const CSRGraph & b) : ia(b.ia), ja(b.ja), tmpfile(b.tmpfile) {}
	CSRGraph & operator=(CSRGraph const & b) { ia = b.ia; ja = b.ja; tmpfile = b.tmpfile; return *this;}
	void LoadBinary(std::ifstream & input)
	{
		size_t header[2];
		input.read(reinterpret_cast<char *>(header),sizeof(size_t)*2);
		ia.resize(header[0]);
		ja.resize(header[1]);
		input.read(reinterpret_cast<char *>(&ia[0]),sizeof(idx_t)*ia.size());
		input.read(reinterpret_cast<char *>(&ja[0]),sizeof(idx_t)*ja.size());
	}
	void LoadBinary(std::string file)
	{
		std::ifstream input(file.c_str(), std::ios::binary);
		LoadBinary(input);
		input.close();
	}
	void SaveBinary(std::ofstream & output) const
	{
		size_t header[2] = {ia.size(),ja.size()};
		output.write(reinterpret_cast<const char *>(header),sizeof(size_t)*2);
		output.write(reinterpret_cast<const char *>(&ia[0]),sizeof(idx_t)*ia.size());
		output.write(reinterpret_cast<const char *>(&ja[0]),sizeof(idx_t)*ja.size());
	}
	void SaveBinary(std::string file) const
	{
		std::ofstream output(file.c_str(), std::ios::binary);
		SaveBinary(output);
		output.close();
	}
	bool Load(std::string file)
	{
		idx_t N,M,L;
		idx_t row, col, last_row;
		double val;
		std::ifstream input(file.c_str());
		if( input.fail() )
		{
			std::cout << "Cannot open " << file << std::endl;
			return false;
		}
		ia.clear();
		ja.clear();
		while( !input.eof() )
		{
			int c = input.peek();
			if( isspace(c) ) continue;
			if( c == '%' ) 
			{
				input.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
				continue;
			}
			input >> N >> M >> L;
			ia.resize(N+1);
			ja.reserve(L);
			break;
		}
		if( ia.empty() ) 
		{
			std::cout << "Error: MatrixMarket header not found." << std::endl;
			return false;
		}
		last_row = 0;
		ia[0] = 0;
		while( !(input >> row >> col >> val).eof() )
		{
			if( !(row == last_row || row == last_row+1) )
			{
				std::cout << "Error: matrix rows are unordered." << std::endl;
				return false;
			}
			if( row != last_row )
			{
				ia[row] = ia[last_row];
				last_row = row;
			}
			ia[row]++;
			ja.push_back(col-1);
		}
		input.close();
		return true;
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
				output << i+1 << " " << ja[j]+1 << " " << 1.0 << std::endl;
		}
		output.close();
	}
	idx_t Columns(bool square = false) const
	{
		idx_t Cols = 0;
		if( square )
			Cols = Size();
		else
		{
			for (idx_t i = 0; i < (idx_t)ja.size(); ++i)
				Cols = std::max(Cols, ja[i] + 1);
		}
		return Cols;
	}
	CSRGraph Transpose(bool square = false) const
	{
		idx_t Cols = Columns(square);
		std::vector< idx_t > iat(Cols+1), ist(Cols,0);
		std::vector< idx_t > jat(Nonzeros()), pos(Nonzeros());
		//find out row size
		for(idx_t i = 0; i < Nonzeros(); ++i)
			pos[i] = ist[ja[i]]++;
		//setup ia array and set writing position
		iat[0] = 0;
		for(idx_t i = 0; i < Cols; ++i)
			iat[i+1] = iat[i] + ist[i];
		for(idx_t i = 0; i < Size(); ++i)
		{
			for(idx_t j = ia[i]; j < ia[i+1]; ++j)
				jat[iat[ja[j]] + pos[j]] = i;
		}
		return CSRGraph(iat,jat);
	}
	// C = A * B
	CSRGraph operator *(const CSRGraph & B) const
	{
		const CSRGraph & A = *this;
		idx_t ColsA = A.Columns(), ColsB = B.Columns();
		//number of columns in A may be less then size of B, i.e. if A is singular
		if(ColsA > B.Size()) throw "Wrong argument size";
		CSRGraph ret;
		RowAccumulatorGraph List(ColsB);
		for(idx_t i = 0; i < A.Size(); ++i)
		{
			for(idx_t j = A.ia[i]; j < A.ia[i+1]; ++j)
				List.Add(B.ia[A.ja[j]],B.ia[A.ja[j]+1],B.ja);
			List.Get(ret.get_ja());
			ret.FinalizeRow();
			List.Clear();
		}
		return ret;
	}
	// C = A + B
	CSRGraph operator +(const CSRGraph & B) const
	{
		const CSRGraph & A = *this;
		//number of columns may be different, i.e. one matrix is singular
		idx_t ColsA = A.Columns(), ColsB = B.Columns(); 
		if(A.Size() != B.Size()) throw "Wrong argument size";
		CSRGraph ret;
		RowAccumulatorGraph List(std::max(ColsA,ColsB));
		for(idx_t i = 0; i < A.Size(); ++i)
		{
			List.Add(A.ia[i],A.ia[i+1],A.ja);
			List.Add(B.ia[i],B.ia[i+1],B.ja);
			List.Get(ret.get_ja());
			ret.FinalizeRow();
			List.Clear();
		}
		return ret;
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
				std::cout << std::setw(w) << "*";
			}
			while(k++ < Cols) std::cout << std::setw(w) << " ";
			std::cout << "|" << std::endl;
		}
	}
	void PrintContents() const
	{
		std::cout << "ia: ";
		for(idx_t i = 0; i < ia.size(); ++i) std::cout << ia[i] << " ";
		std::cout << std::endl;
		std::cout << "ja: ";
		for(idx_t i = 0; i < ja.size(); ++i) std::cout << ja[i] << " ";
		std::cout << std::endl;
	}
	
	void SortColumns()
	{
		for(idx_t k = 0; k < Size(); ++k)
		{
			idx_t i, j, n = ia[k+1] - ia[k], s = ia[k];  
			for (i = 0; i < n-1; i++)
			{
				for (j = 0; j < n-i-1; j++)
				{
					if (ja[s+j] > ja[s+j+1])
						std::swap(ja[s+j],ja[s+j+1]);
				}
			}
		}
	}
	bool Symmetric() const
	{
		idx_t normm = 0;
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
					normm++;
				else if(i1 != i)
				{
					for (idx_t j1 = ia[i1]; j1 < ia[i1+1]; j1++) if (ja[j1] == i)
						mrk[j1] = !mrk[j1];
				}
			}
		}
		
		if( normm == 0 )
			return true;
		return false;
	}
	bool Compare(CSRGraph const & b) const
	{
		bool ret = true;
		if( Size() == b.Size() )
		{
			for(idx_t i = 0; i < Size(); ++i)
			{
				if( RowSize(i) != b.RowSize(i) )
				{
					std::cout << "At row "<<i<< " row sizes are " << RowSize(i) << " and " << b.RowSize(i) << std::endl;
					ret = false;
				}
				else
				{
					for(idx_t j = ia[i]; j < ia[i+1]; ++j)
					{
						ret = false;
						for(idx_t j1 = b.ia[i]; j1 < b.ia[i+1]; ++j1) if( ja[j] == b.ja[j1] )
						{
							ret = true;
							break;	
						}
						if( !ret ) 
							std::cout << "At ("<<i<<","<<ja[j]<<") no matching column in other graph" << std::endl;
					}
				}
			}
			if( ret ) 
				std::cout << "Graphs are similar" << std::endl;
			else
				std::cout << "Graphs are different" << std::endl;
		}
		else 
		{
			std::cout << "Graph sizes " << Size() << " and " << b.Size() << std::endl;
			ret = false;
		}
		return ret;
	}
	size_t Nonzeros() const {return ja.size();}
	idx_t Size() const {return ia.size()-1;}
	idx_t RowSize(idx_t i) const {return ia[i+1]-ia[i];}
	idx_t Col(idx_t i, idx_t k) const {return ja[ia[i]+k];}
	idx_t Col(const std::pair<idx_t,idx_t> & p) const {return Col(p.first,p.second);}
	idx_t & Col(idx_t i, idx_t k) {return ja[ia[i]+k];}
	idx_t & Col(const std::pair<idx_t,idx_t> & p) {return Col(p.first,p.second);}
	void ConcatRow(const std::vector<idx_t> & ja_part)
	{
		Insert(ja_part);
		FinalizeRow();
	}
	void Insert(const std::vector<idx_t> & ja_part) {ja.insert(ja.end(),ja_part.begin(),ja_part.end());}
	void PushBack(idx_t col) { ja.push_back(col); }
	void FinalizeRow() { ia.push_back(ja.size()); }
	void Clear() 
	{
		if( Offloaded() )
		{
			std::remove(tmpfile.c_str());
			tmpfile = "";
		}
		ia.clear(); 
		ja.clear(); 
		ia.push_back(0);
	}
	size_t Bytes() const {return get_bytes(ia) + get_bytes(ja);}
	bool Offloaded() const
	{
		if( tmpfile != "" )
			return true;
		return false;
	}
	void Offload()
	{
		if( Size() != 0 )
		{
			//get temporary file name
			tmpfile = std::tmpnam(NULL);
			std::cout << __FILE__ << ":" << __LINE__ << " Info: offloading graph to " << tmpfile << "." << std::endl;
			SaveBinary(tmpfile);
			Clear();
			Shrink();
		}
		else std::cout << __FILE__ << ":" << __LINE__ << " Warning: requested offloading for empty graph." << std::endl;
	}
	void Reload()
	{
		if( Size() == 0 )
		{
			if( tmpfile != "" )
			{
				LoadBinary(tmpfile);
				std::remove(tmpfile.c_str());
				tmpfile = "";
			}
			else
				std::cout << __FILE__ << ":" << __LINE__ << " Error: no temporary file with graph data." << std::endl;
		} else std::cout << __FILE__ << ":" << __LINE__ << " Error: reloading data to non-empty graph." << std::endl;
	}
	void ReserveSize(size_t size) {ia.reserve(size+1);}
	void ReserveNonzeros(size_t nonzeros) {ja.reserve(nonzeros);}
	void Shrink() { ia.shrink_to_fit(); ja.shrink_to_fit(); }
	
	std::vector<idx_t> & get_ia() {return ia;}
	std::vector<idx_t> & get_ja() {return ja;}
	const std::vector<idx_t> & get_ia() const {return ia;}
	const std::vector<idx_t> & get_ja() const {return ja;}
};

#endif //_GRAPH_H
