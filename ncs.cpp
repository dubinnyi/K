#include<vector>
#include<set>
#include<map>
#include<string>

using namespace std; 

class labeltype; 

class spectrum {
public:
	string name;
	spectrum(string sname) {
		name = sname;
	}

	int has_signal(labeltype label_type_1, labeltype label_type_2) {
		if (name == "HSQC") {
			return ((label_type_2.isotopes / 100) % 10);
		}
		else if(name == "HNCO"){
			return ((label_type_2.isotopes / 100) % 10)*(label_type_1.isotopes % 10);
		}
		else if (name == "HNCA") {
			return ((label_type_2.isotopes / 100) % 10)*(((label_type_2.isotopes / 10) % 10) + ((label_type_1.isotopes / 10) % 10));
		}
		else if (name == "HNCOCA") {
			return ((label_type_2.isotopes / 100) % 10)*((label_type_1.isotopes / 10) % 10)*(label_type_1.isotopes % 10);
		}
		else if (name == "COfHNCA") {
			int k;
			if (label_type_1.isotopes % 10 == 0) 
				k = 1;
			else
			k = 0;

			return k * ((label_type_1.isotopes / 10) % 10)*((label_type_2.isotopes / 100) % 10);
		}
		else if (name == "DQHNCA") {
			return ((label_type_2.isotopes / 100) % 10)*(((label_type_2.isotopes / 10) % 10) * ((label_type_1.isotopes / 10) % 10));
		}

		else if (name == "HNCACO") {
			return ((label_type_2.isotopes / 100) % 10)*(((label_type_2.isotopes / 10) % 10) * ((label_type_2.isotopes% 10)));

		}

		else if (name == "HNCAfCO") {
			int k;
			if ((label_type_2.isotopes / 10)%10 == 0)
				k = 1;
			else
				k = 0;

			return ((label_type_2.isotopes / 100) % 10)*k * ((label_type_2.isotopes % 10));
		}

		

	}


};

class labeltype {
public:
	string name;
	int  isotopes;
	labeltype(string lname,int  lisotopes) {
		name = lname;
		isotopes = lisotopes; 
	}
};

class constants {
	constants() {
		spectrum HSQC("HSQC");
		spectrum HNCO("HNCO");
		spectrum HNCA("HNCA");
		spectrum HNCOCA("HNCOCA");
		spectrum DQHNCA("DQHNCA");
		spectrum COfHNCA("COfHNCA");
		spectrum HNCACO("HNCACO");
		spectrum HNCAfCO("HNCAfCO");
		const vector <spectrum> basic_spectra = { HSQC, HNCO, HNCA, HNCOCA, DQHNCA, COfHNCA, HNCACO, HNCAfCO };
		labeltype typeX("X", 000);
		labeltype typeN("N", 100);
		labeltype typeC("C", 001);
		labeltype typeD("D", 111);
		labeltype typeA("A", 010);
		labeltype typeT("T", 110);
		labeltype typeS("S", 101);
		labeltype typeF("F", 011);
		const vector <labeltype> BASIC_TYPES = {typeX, typeN, typeC, typeD, typeA, typeT, typeS, typeF	};
	}
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

	//name, spectra_list, label_types, deuterated = False)
	NCS(string name_ncs, vector<spectrum>spectra_list_ncs, vector<labeltype> label_types_ncs, bool deuterated_ncs = 0) {
		name = name_ncs;
		spec_list = spectra_list_ncs;
		label_types = label_types_ncs;
		deuterated = deuterated_ncs;
		for (labeltype l : label_types) {
			label_dict[l.name] = l;
		}
		letters = {"a", "b", "c", "d", "e", "f", "g", "h", "i", "j"};
		for (spectrum s : spec_list) {
			//vectors.push_back([0]);
			vectors.push_back({0});
		}
	}

	void make_coding_table() {
		vector <vector <string>> codes_table;
		for (int j = 0; j < label_types.size(); j++) {
			for (int i = 0; i < label_types.size(); i++) {
				codes_table.push_back({ "0" });
			}
		}
		vector<int> vec;
		vector<string> result;
		set<string> s;
		int power;
		int code;
		//labeltype label_2, label_1;
		for (int i = 0; i < label_types.size(); i++) {
			labeltype label_2 = label_types[i];

			for (int j = 0; j < label_types.size(); j++) {
				labeltype label_1 = label_types[j];
				vec = {};
				for (spectrum spect : spec_list) {
					vec.push_back(spect.has_signal(label_1, label_2));
				}
				if (find(vectors.begin(), vectors.end(), vec) != vectors.end()) {
					code = distance(vectors.begin(), find(vectors.begin(), vectors.end(), vec));
				}
				else {
					code = vectors.size();
					vectors.push_back(vec);
				}
				if (code > 9) {
					codes_table[j][i] = letters[code - 10];
				}
				else {
					codes_table[j][i] = to_string(code);
				}
			}
		}
		for (int i = 0; i < label_types.size(); i++) {
			result = {};
			for (vector <string> row : codes_table) {
				result.push_back(row[i]);
			}
			set <string> s(result.begin(), result.end());
			power = s.size();
			label_power[label_types[i]] = power;
			if (power == 1 and find(NITRO_TYPES.begin(), NITRO_TYPES.end(), label_types[i]) != NITRO_TYPES.end()) {

				throw ("Error");
			}
		}
		map <labeltype, string> subdict;
		for (int i = 0; i < label_types.size(); i++) {
			labeltype  label_1 = label_types[i];

			for (int j = 0; j < label_types.size(); j++) {
				labeltype label_2 = label_types[j];
				subdict[label_2] = codes_table[i][j];
		}
			codes_dict[label_1] = subdict;
	}

	}



};

#endif // NCS_H_INCLUDED
