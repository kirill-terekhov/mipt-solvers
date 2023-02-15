#ifndef _VECTOR_H
#define _VECTOR_H
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <cmath>

void LoadVectorBinary(std::string file, std::vector<double> & v)
{
	std::ifstream input(file.c_str(),std::ios::binary);
	size_t header[1];
	input.read(reinterpret_cast<char *>(header),sizeof(size_t));
	v.resize(header[0]);
	input.read(reinterpret_cast<char *>(&v[0]),sizeof(double)*v.size());
	input.close();
}

void SaveVectorBinary(std::string file, const std::vector<double> & v)
{
	std::ofstream output(file.c_str(),std::ios::binary);
	const size_t header[1] = {v.size()};
	output.write(reinterpret_cast<const char *>(header),sizeof(size_t));
	output.write(reinterpret_cast<const char *>(&v[0]),sizeof(double)*v.size());
	output.close();
}

void LoadVector(std::string file, std::vector<double> & v)
{
	idx_t N, row;
	double val;
	std::ifstream input(file.c_str());
	if( input.fail() )
	{
		std::cout << "Cannot open " << file << std::endl;
		exit(-1);
	}
	v.clear();
	while( !input.eof() )
	{
		int c = input.peek();
		if( isspace(c) ) continue;
		if( c == '%' )
		{
			input.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
			continue;
		}
		input >> N;
		v.resize(N);
		break;
	}
	if( v.empty() ) 
		throw "Vector is empty";
	row = 0;
	while( !(input >> val).eof() )
		v[row++] = val;
	input.close();
}

void SaveVector(std::string file, const std::vector<double> & v)
{
	std::ofstream output(file.c_str());
	output << v.size() << std::endl;
	output << std::scientific;
       	output.precision(16);
	for(idx_t i = 0; i < (idx_t)v.size(); ++i)
		output << v[i] << std::endl;
	output.close();
}

void SaveVector(std::string file, const std::vector<idx_t> & v)
{
	std::ofstream output(file.c_str());
	output << v.size() << std::endl;
	for(idx_t i = 0; i < (idx_t)v.size(); ++i)
		output << v[i] << std::endl;
	output.close();
}

template<typename vtype>
void SaveVector(std::string file, const std::vector<vtype>& v)
{
	std::ofstream output(file.c_str());
	output << v.size() << std::endl;
	for (idx_t i = 0; i < (idx_t)v.size(); ++i)
		output << v[i] << std::endl;
	output.close();
}


INLINE double Dot(const std::vector<double> & x, const std::vector<double> & y)
{
	if( x.size() != y.size() ) throw "Wrong argument size";
	double norm = 0;
	for(idx_t i = 0; i < (idx_t)x.size(); ++i)
		norm += x[i]*y[i];
	return norm;
}


INLINE double Norm(const std::vector<double> & x)
{
	return sqrt(Dot(x,x));
}

void Zero(std::vector<double> & x)
{
	std::fill(x.begin(),x.end(),0.0);
}

#endif //_VECTOR_H
