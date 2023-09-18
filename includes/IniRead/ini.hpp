#ifndef INI_H__
#define INI_H__

#include<map>
#include<string>
#include<sstream>
#include<fstream>

template<class T> bool ini_read(const char * fname, T & value, const char * label, const T & default_=T() ){
	std::map<std::string, std::string> m;
	std::ifstream ifs(fname);
	std::string label_, value_;
	while(ifs >> label_ && getline(ifs, value_) ) 
		m.insert(std::make_pair(label_, value_));
			
	std::map<std::string, std::string>::iterator p=m.find(label);
	if ( p!=m.end() ){
		std::istringstream os(m[label]);
		os >> value;
		return true;
	}
	value=default_;
	return false;
}

template<class T> T ini_read(const char * fname, const char * label, const T & default_=T()){
	T tmp;
	ini_read(fname, tmp, label, default_);
	return tmp;
}

#endif 

