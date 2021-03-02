#ifndef _PARAMS_H
#define _PARAMS_H

#define NOMINMAX
#include <sstream>
#include <map>

class Parameters
{
	std::map<std::string, std::string> list;
	std::map<std::string, Parameters> sub_list;
	template<typename RetType>
	static RetType from_string(std::string value)
	{
		RetType ret;
		value = value.substr(0,value.find("#")); //erase comment
		std::istringstream convert(value.c_str());
		convert >> ret;
		return ret;
	}
	template<typename ValType>
	static std::string to_string(ValType value)
	{
		std::stringstream convert;
		convert << value;
		return convert.str();
	}
	static std::string strip_spaces(std::string value)
	{
		std::string::const_iterator it = value.begin();
		std::string::const_reverse_iterator rit = value.rbegin();
		while (it != value.end() && isspace(*it)) it++;
		while (rit.base() != it && isspace(*rit)) rit++;
		return std::string(it, rit.base());
	}
	static std::string get_prefix(std::string name)
	{
		std::size_t pos = name.find(":");
		if( pos != std::string::npos )
			return name.substr(0,pos);
		return "";
	}
	static std::string strip_prefix(std::string name)
	{
		std::size_t pos = name.find(":");
		if( pos != std::string::npos )
			return name.substr(pos+1);
		return name;
	}
	static std::ostream & write_tabs(std::ostream & stream, int tabs)
	{
		for(int k = 0; k < tabs; ++k) stream << "\t";
		return stream;
	}
	void pretty_print(std::ostream & output, int tabs) const
	{
		for(std::map<std::string,std::string>::const_iterator it = list.begin(); it != list.end(); ++it)
			write_tabs(output,tabs) << it->first << " = " << it->second  << std::endl;
		for(std::map<std::string,Parameters>::const_iterator it = sub_list.begin(); it != sub_list.end(); ++it)
		{
			write_tabs(output,tabs) << it->first << ":" << std::endl;
			it->second.pretty_print(output,tabs+1);
			write_tabs(output,tabs) << "/" << std::endl;
		}
	}
	std::istream & get_name(std::istream & input, std::string & ret) const
	{
		char c;
		ret.clear();
		while( input.get(c) && c != '\n' && c != '=' )
			ret.push_back(c);
		ret = strip_spaces(ret);
		if( c == '=' )
			input.unget();
		return input;
	}
	void pretty_read(std::istream & input)
	{
		char equal;
		std::string name,value;
		while( !get_name(input,name).eof() )
		{
			//~ std::cout << name;
			if (name.empty()) continue;
			if( name[name.size()-1] == ':' )
				sub_list[name.substr(0,name.size()-1)].pretty_read(input);
			else if( name == "/" )
				break;
			else if( name != "" )
			{
				while( input.get(equal) && equal != '=' );
				if( equal != '=' ) throw "Unexpected character";
				//~ input >> value;
				if( std::getline(input,value,'\n') )
					Set(name,value);
				//~ input.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
			}
				
		}
	}
public:
	Parameters() {}
	Parameters(const Parameters & b) : list(b.list), sub_list(b.sub_list) {}
	Parameters & operator = (Parameters const & b) {list = b.list; sub_list = b.sub_list; return *this;}
	Parameters & SubParameters(std::string prefix) {return sub_list[prefix];}
	const Parameters & SubParameters(std::string prefix) const
	{
		std::map<std::string,Parameters>::const_iterator find = sub_list.find(prefix);
		if( find == sub_list.end() )
			throw "Invalid prefix name";
		return find->second;
	}
	bool SubParametersExist(std::string prefix) const
	{
		if( sub_list.find(prefix) == sub_list.end() )
			return false;
		return true;
	}
	template<typename ValType>
	Parameters & SubParametersSearchRecursive(std::string prefix, std::string param, const ValType & value)
	{
		std::string cmp = to_string(value);
		Parameters * p = this;
		while( p->SubParametersExist(prefix) )
		{
			if( p->Have(param) && p->Get<std::string>(param) == cmp )
				break;
			p = &p->SubParameters(prefix);
		}
		return *p;
	}
	template<typename ValType>
	const Parameters & SubParametersSearchRecursive(std::string prefix, std::string param, const ValType & value) const
	{
		std::string cmp = to_string(value);
		const Parameters * p = this;
		while( p->SubParametersExist(prefix) )
		{
			if( p->Have(param) && p->Get<std::string>(param) == cmp )
				break;
			p = &p->SubParameters(prefix);
		}
		return *p;
	}
	//~ void ConnectParameters(std::string prefix, const Parameters & params) { sub_list[prefix] = params; }
	void Print(std::ostream & stream = std::cout, std::string prefix = "") const
	{
		for(std::map<std::string,std::string>::const_iterator it = list.begin(); it != list.end(); ++it)
		{
			if( prefix != "" ) stream << prefix << ":";
			stream << it->first << " = " << it->second << std::endl;
		}
		for(std::map<std::string,Parameters>::const_iterator it = sub_list.begin(); it != sub_list.end(); ++it)
		{
			if( prefix != "" )
				it->second.Print(stream,prefix+":"+it->first);
			else
				it->second.Print(stream,it->first);
		}
	}
	void PrettyPrint() const
	{
		std::cout << "Method:" << std::endl;
		pretty_print(std::cout,1);
		std::cout << "/" << std::endl;
	}
	void SaveRaw(std::string file) const
	{
		std::ofstream output(file.c_str());
		if( output.fail() )
		{
			std::cout << "Cannot open " << file << std::endl;
			return;
		}
		Print(output);
		output.close();
	}
	void LoadRaw(std::string file)
	{
		std::string name, value, prefix;
		char equal;
		std::ifstream input(file.c_str());
		if( input.fail() )
		{
			std::cout << "Cannot open " << file << std::endl;
			return;
		}
		while(!(input >> name >> equal >> value).eof())
		{
			if( equal != '=' ) throw "Unexpected character";
			if( value != "*" ) Set(name,value);
			input.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
		}
		input.close();
	}
	void Save(std::string file)
	{
		std::ofstream output(file.c_str());
		if( output.fail() )
		{
			std::cout << "Cannot open " << file << std::endl;
			return;
		}
		output << "Method:" << std::endl;
		pretty_print(output,1);
		output << "/" << std::endl;
		output.close();
	}
	void Load(std::string file)
	{
		std::string method;
		std::ifstream input(file.c_str());
		if( input.fail() )
		{
			std::cout << "Cannot open " << file << std::endl;
			return;
		}
		input >> method;
		if( method == "Method:" )
			pretty_read(input);
		else std::cout << "Parameters file should start with \"Method:\"" << std::endl;
		input.close();
	}
	template<typename ValType> void SetRecursive(std::string name, const ValType & value) 
	{
		if( Have(name) )
			list[name] = to_string(value); 
		for(std::map<std::string,Parameters>::iterator it = sub_list.begin(); it != sub_list.end(); ++it)
			it->second.SetRecursive(name,value);
	}
	template<typename ValType> void Set(std::string name, const ValType & value) 
	{
		std::string prefix;
		Parameters * p = this;
		while( (prefix = get_prefix(name)) != "" )
		{
			p = &p->SubParameters(prefix);
			name = strip_prefix(name);
		}
		p->list[name] = strip_spaces(to_string(value));
	}
	template<typename RetType> RetType Get(std::string name) const 
	{
		std::string prefix;
		const Parameters * p = this;
		while( (prefix = get_prefix(name)) != "" )
		{
			p = &p->SubParameters(prefix);
			name = strip_prefix(name);
		}
		std::map<std::string,std::string>::const_iterator find = p->list.find(name);
		if( find == p->list.end() )
			throw "Invalid parameter name";
		return from_string<RetType>(find->second); 
	}
	bool Have(std::string name) const {return list.find(name) != list.end();}
};

#endif //_PARAMS_H
