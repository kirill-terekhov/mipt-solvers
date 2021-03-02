#ifndef _METHOD_H
#define _METHOD_H
#include "csrmatrix.h"
#include "params.h"
#include <sstream>


class Methods
{
	Parameters params;
public:
	static Parameters DefaultParameters() {return Parameters();}
	Parameters & GetParameters() {return params;}
	const Parameters & GetParameters() const {return params;}
	void SetParameters(const Parameters & p) {params = p;}
	virtual bool Setup(const CSRMatrix & A) = 0;
	virtual bool Solve(const std::vector<double> & b, std::vector<double> & x) const = 0;
	virtual size_t Bytes() const = 0;
	virtual ~Methods() {}
};

class DummySolver : public Methods
{
public:
	DummySolver() {GetParameters() = DefaultParameters();}
	static Parameters DefaultParameters() 
	{
		Parameters ret;
		ret.Set("name","DummySolver");
		return ret;
	}
	bool Setup(const CSRMatrix & A) {(void)A; return true;}
	bool Solve(const std::vector<double> & b, std::vector<double> & x) const {x = b; return true;}
	size_t Bytes() const {return 0;}
};

template<typename Solver>
class DummyPreprocessor : public Methods
{
	Solver S;
public:
	static Parameters DefaultParameters() 
	{
		Parameters ret; 
		ret.Set("name","DummyPreprocessor");
		ret.Set("level","*");
		ret.SubParameters("Solver") = Solver::DefaultParameters(); 
		return ret;
	}
	DummyPreprocessor() : S() {GetParameters() = DefaultParameters();}
	bool Setup(const CSRMatrix & A) {S.SetParameters(GetParameters().SubParameters("Solver")); return S.Setup(A);}
	bool Solve(const std::vector<double> & b, std::vector<double> & x) const {return S.Solve(b,x);}
	size_t Bytes() const {return S.Bytes();}
};

template<typename T>
static std::string to_string(const T & input)
{
	std::ostringstream stream;
	stream << input;
	return stream.str();
}

#endif //_METHOD_H
