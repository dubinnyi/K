#ifndef NCS_H_INCLUDED
#define NCS_H_INCLUDED

#include<string>
#include<vector>
#include<map>

using namespace std; 

class labeltype {
public:
	string name;
	int  isotopes;
  labeltype();
  labeltype(const labeltype& t);
	labeltype(string lname, int  lisotopes);
//	auto operator<=>(const labeltype&) = default
	bool operator<(const labeltype & t2);
};
	
bool operator<(const labeltype& t1, const labeltype& t2);
bool operator==(const labeltype& t1, const string& s2);
bool operator==(const string& s1, const labeltype& t2);


class spectrum {
public:
	string name;
	spectrum(string sname);
	int has_signal(labeltype label_type_1, labeltype label_type_2);
};

class NCS {
public:
	const vector <string> NITRO_TYPES = {"N", "D", "S", "T"};
	string name;
	vector<spectrum>spec_list;
	vector<labeltype> label_types;
	bool deuterated;
	map<string, labeltype> label_dict;
	vector <string> letters; 
	//vector <> spectra_numbers; 
	map<labeltype, int> label_power; 
	map <labeltype, map <labeltype, string>> codes_dict;
	vector <vector <int>> vectors;
	map <labeltype, string> subdict;

	//name, spectra_list, label_types, deuterated = False)
	NCS(string name_ncs, vector<spectrum>spectra_list_ncs, vector<labeltype> label_types_ncs, bool deuterated_ncs = false);
  void make_coding_table(void);
};

#endif // NCS_H_INCLUDED
