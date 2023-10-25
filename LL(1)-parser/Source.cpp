#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <algorithm>
using namespace std;
class Grammar {
	vector <string> terminals;
	vector <string> non_terminals;
	string axiom;
	map<string, map<int, vector<string>>> rules; // map<left side of the rule(non-terminal), map<number of the rule for the non-terminal from the left side of the rule, vector<string> of terminals/non-terminals for this rule>>
	/*
	Example: A: + B A | eps
	map<string, map<int, vector<string>>> let string=S, int=I, vector<string>=Strings
			   => I=0 && Strings={"+", "B", "A"}
			  |	
	Then S="A" 
			  |
			   => I=1 && Strings={"eps"}
	*/
	//map<string, set<string>> first_k; //map for the first_k sets; OLD FUNCTIONALITY
	map <string, map<string, vector<int>>> first_k2;
	/*
	map<string, map<string,vector<int>>> where string is the name of non-terminal; 
	in map<string, vector<int>> string is a string of terminals of length<=k which can be derived from non-terminal, 
	vector<int> is a vector of positions of the ends of the terminal characters int the string

	map <string1, map<string2, vector<int>>> let string1=S1, int=I, vector<int>=Ints, string2=S2
	Example:
	
									   =>S2="id*id" && Ints={2,1,2}, because len("id")=2 and len("*")=1
									  |  
									  |  
	  then at some step assume S1="A"  =>S2="(*id" && Ints={1,1,2}, because len("(")=len("*")=1 && len("id")=2 and "(","*","id" are different terminals
									  |
									  |
									   =>S2="id)" && Ints={2,1}, because len("id")=2 and len(")")=1

	With this structure we can add (in the sense of a binary operation on dictionary sets) terminal symbols whose length is >1
	*/
	set<string> epsilon; //for epsilon non-terminals
	set<string> leftRecursive; //for left-recursive non-terminals 

public:

	void read(string filename) {
		//reading axiom
		string axiom;
		string delimiter_axiom = "%start ";
		ifstream file(filename);
		getline(file, axiom);
		axiom.erase(0, axiom.find(delimiter_axiom) + delimiter_axiom.length());
		this->axiom = axiom;

		//reading tokens
		string token;
		string delimiter_token = "%token ";
		int pos;
		getline(file, token);
		token.erase(0, token.find(delimiter_token) + delimiter_token.length());
		while ((pos = token.find(" ")) != string::npos) {
			this->terminals.push_back(token.substr(0, pos));
			token.erase(0, pos + 1);
		}
		this->terminals.push_back(token);

		//reading rules
		while (!file.eof()) {
			string rule,non_terminal;
			string delimiter_transform = ":", delimiter_sym = " ", delimiter_or = "|";
			vector <string> symbols; //for terminals + non-terminals of one rule
			map <int, vector<string>> rules; // for the rules of one non-terminal
			int rules_count = 0;
			int pos_or;
			getline(file, rule); // rules line such as  A: + B A | eps
			pos_or = rule.find(delimiter_or);
			non_terminal = rule.substr(0, rule.find(delimiter_transform)); //left side of the rule
			this->non_terminals.push_back(non_terminal);
			rule.erase(0, rule.find(delimiter_sym)+1);
			if (rule.find(delimiter_sym) == string::npos) { symbols.push_back(rule); rules[rules_count] = symbols; this->rules[non_terminal] = rules; } //if there's only one terminal/non-terminal
			else {
				while (rule.find(delimiter_sym) != string::npos) { //while we have terminals/non-terminals
					while (pos_or != 0 && rule.find(delimiter_sym) != string::npos) { //while we have rules for our non-terminal from the left side
						pos = rule.find(delimiter_sym);
						symbols.push_back(rule.substr(0, pos)); // add terminal/non-terminal
						rule.erase(0, pos + 1); //cut the string
						pos_or = rule.find(delimiter_or);
					}
					if (rule.find(delimiter_sym) == string::npos) symbols.push_back(rule); //if there's the last terminal/non-terminal 
					rules[rules_count] = symbols;
					++rules_count;
					symbols.erase(symbols.begin(), symbols.end());
					pos_or = rule.find(delimiter_or);
					if (pos_or != string::npos) { //if we still have rules to parse
						rule.erase(0, rule.find(delimiter_sym) + 1); //cut the string
						pos_or = rule.find(delimiter_or);
						if (rule.find(delimiter_sym) == string::npos) { symbols.push_back(rule); rules[rules_count] = symbols; symbols.erase(symbols.begin(), symbols.end()); }
					}
					
				}
				this->rules[non_terminal] = rules;
			}
		}
		file.close();
	}

	bool isNonTerminal(string s) {
		for (int i = 0; i < non_terminals.size(); ++i) if (s == non_terminals[i]) return true;
		return false;
	}

	bool isTerminal(string s) {
		for (int i = 0; i < terminals.size(); ++i) if (s == terminals[i]) return true;
		return false;
	}

	//OLD FUNCTIONALITY
	//sum binary operation on the words of a language, k is the length
	/*set<string> sumOfSets(set<string> set1, set<string> set2, int k) {
		set<string> res;
		for (string x : set1) {
			for (string y : set2) {
				if (x == "eps") res.insert(y.substr(0, k));
				else if (y == "eps") res.insert(x.substr(0, k));
				else res.insert((x + y).substr(0, k));
			}
		}
		return res;
	}*/

	// first_k building, size is the length
	/*void build_first_k(int size) {
		set<string> temp;
		string temp_term;
		vector<string> rule; //vector for the rule (one of the right parts)
		bool all_is_terminals;
		for (int i = 0; i < non_terminals.size(); ++i) {
			first_k[non_terminals[i]] = temp; //creating sets for each non-terminal
		}
		for (int i = 0; i < non_terminals.size(); ++i) {
			int k = rules[non_terminals[i]].size(); // number of rules for the non-terminal
			for (int j = 0; j < rules[non_terminals[i]].size();++j) { //we check each rule
				rule = rules[non_terminals[i]][j];
				all_is_terminals = true;
				for (int b = 0; b < rule.size(); ++b) {
					if (rule[b] != "eps") {
						if (isNonTerminal(rule[b])) {
							all_is_terminals = false;
							break;
						}
					}
				}
				//if the rule consist only of terminals then we add it the the first_k(non-terminal[i])
				if (all_is_terminals) {
					temp_term = "";
					for (int count = 0; count < size; ++count) { if (count >= rule.size()) break; temp_term += rule[count]; }
					first_k[non_terminals[i]].insert(temp_term);
				}
			}
		}
		map<string, set<string>> temp_configuration=first_k; //to check when the iterative algorithm becomes stable
		bool configuration = false,set_is_empty;
		set<string> terminal_symbol, set_for_nonterminal;
		while (configuration == false) { //while the algorithm isn't stable
			for (int i = 0; i < non_terminals.size(); ++i) { //for each non-terminal
				int k = rules[non_terminals[i]].size(); //number of rules for the non-terminal
				for (int j = 0; j < rules[non_terminals[i]].size(); ++j) { //for each rule
					rule = rules[non_terminals[i]][j];
					set_for_nonterminal.clear(); //for iterative union of sets
					if (isTerminal(rule[0])) set_for_nonterminal.insert(rule[0]); //if the first sym is terminal we add it to the set
					else set_for_nonterminal = first_k[rule[0]]; //otherwise we add a set first_k(non-terminal[i])
					if (set_for_nonterminal.empty()) continue;
					else {
						set_is_empty = false;
						for (int b = 1; b < rule.size(); ++b) { //for each terminal/non-terminal in the rule
							//getting previous first_k(non-terminal[i]) set or terminal symbol since first_k(terminal)={terminal}
							if (isTerminal(rule[b])) {
								terminal_symbol.clear(); terminal_symbol.insert(rule[b]);
							}
							else terminal_symbol = first_k[rule[b]];
							if (terminal_symbol.empty()) {
								set_is_empty = true;
								break;
							}
							set_for_nonterminal = sumOfSets(set_for_nonterminal, terminal_symbol, size); //iterative union of sets
						}
						if (!set_is_empty) {
							for (string x : set_for_nonterminal) {
								first_k[non_terminals[i]].insert(x);
							}
						}  //if both of the sets weren't empty we unite set on the previous step and on the current
					}
				}

			}
			if (first_k == temp_configuration) configuration = true; // we check whether the algorithm has stabilized
			else temp_configuration = first_k;
		}
		first_k.erase("eps");
	}*/

	//sum binary operation on the words of a language, k is the length
	map<string, vector<int>> sumOfSets(map<string, vector<int>> set1, map<string, vector<int>> set2, int k) {
		map<string, vector<int>>res; //for the result
		vector<int> res3; //for vectors which will be included to the map
 		int iterator = 0;
		string concat; // string of terminal characters
		int size,sum1=0,sum2=0,it,position;
				for (auto const& x : set1) {
					for (auto const& y : set2) {
						if (y.first == "eps" && x.first == "eps") { //both are eps, then x+y=eps
							res3.clear();
							res3.push_back(3);
							concat = "eps";
							res[concat] = res3;
						}
						else if (x.first == "eps") { //x=eps, then x+y=y cutted to k characters
							sum2 = 0;
							res3.clear();
							concat = "";
							for (it = 0; it < k; ++it) {
								if (it < y.second.size()) { sum2 += y.second[it]; res3.push_back(y.second[it]); } //taking <=k terminals to the answer from the string y
							}
							concat = y.first.substr(0, sum2);
							res[concat] = res3;
						}
						else if (y.first == "eps") {  //y=eps, then x+y=x cutted to k characters
							sum1 = 0;
							res3.clear();
							concat = "";
							for (it = 0; it < k; ++it) {
								if (it < x.second.size()) { sum1 += x.second[it]; res3.push_back(x.second[it]); }  //taking <=k terminals to the answer from the string x
							}
							concat = x.first.substr(0, sum1);
							res[concat] = res3;
						}
						else { //both arent equal eps
							iterator = 0;
							concat = "";
							sum1 = sum2 = 0;
							it = 0;
							position = 0;
							res3.clear();
							for (it = 0; it < k; ++it) {
								if (position < x.second.size()) { sum1 += x.second[position]; res3.push_back(x.second[position]); ++position; } //taking <=k terminals to the answer from the string x
								else break;
							}
							position = 0;
							for (it ; it < k; ++it) {
								if (position < y.second.size()) { sum2 += y.second[position]; res3.push_back(y.second[position]); ++position; } // if the length of the result string is still <k taking <=k terminals to the answer from the string y
								else break;
							}
							concat = x.first.substr(0, sum1);
							concat += y.first.substr(0, sum2); //n terminal from the string x and maybe another m terminals from the string y
							res[concat] = res3;
						}
					}
				}
				return res;
	}

	void build_first_k(int size) {
		map<string, vector<int>> temp;
		string temp_term;
		vector<string> rule; //vector for the rule (one of the right parts)
		map<string, vector<int>> map_rules;
		vector<int> prepositions;
		bool all_is_terminals;
		for (int i = 0; i < non_terminals.size(); ++i) {
			first_k2[non_terminals[i]] = temp; //creating maps for each non-terminal
		}
		for (int i = 0; i < non_terminals.size(); ++i) {
			int k = rules[non_terminals[i]].size(); // number of rules for the non-terminal
			for (int j = 0; j < rules[non_terminals[i]].size(); ++j) { //we check each rule
				rule = rules[non_terminals[i]][j];
				all_is_terminals = true;
				for (int b = 0; b < rule.size(); ++b) {
					if (rule[b] != "eps") {
						if (isNonTerminal(rule[b])) {
							all_is_terminals = false;
							break;
						}
					}
				}
				//if the rule consist only of terminals then we add it the the first_k(non-terminal[i])
				if (all_is_terminals) {
					temp_term = "";
					for (int count = 0; count < size; ++count) { if (count >= rule.size()) break; temp_term += rule[count]; }
					map_rules.clear(); prepositions.clear();
					prepositions.push_back(temp_term.length()); //inserting the length of the terminal
					map_rules[temp_term] = prepositions;
					first_k2[non_terminals[i]]=map_rules;
				}
			}
		}
		map<string,map<string, vector<int>>> temp_configuration=first_k2; //to check when the iterative algorithm becomes stable
		bool configuration = false, set_is_empty;
		map<string, vector<int>> terminal_symbol, set_for_nonterminal;
		
		vector<int> positions;
		while (configuration == false) { //while the algorithm isn't stable
			for (int i = 0; i < non_terminals.size(); ++i) { //for each non-terminal
				int k = rules[non_terminals[i]].size(); //number of rules for the non-terminal
				for (int j = 0; j < rules[non_terminals[i]].size(); ++j) { //for each rule
					rule = rules[non_terminals[i]][j];
					set_for_nonterminal.clear(); //for iterative union of sets 
					if (isTerminal(rule[0])) { positions.clear(); positions.push_back(rule[0].length());set_for_nonterminal[rule[0]]=positions;   }//if the first sym is terminal we add it to the map
					else set_for_nonterminal = first_k2[rule[0]]; //otherwise we add a set first_k(non-terminal[i])
					if (set_for_nonterminal.empty()) continue;
					else {
						set_is_empty = false;
						for (int b = 1; b < rule.size(); ++b) { //for each terminal/non-terminal in the rule
							//getting previous first_k(non-terminal[i]) set or terminal symbol since first_k(terminal)={terminal}
							if (isTerminal(rule[b])) {
								terminal_symbol.clear();
								positions.clear();
								positions.push_back(rule[b].length()); 
								terminal_symbol[rule[b]]=positions;
								positions.clear(); 
							}
							else terminal_symbol = first_k2[rule[b]];
							if (terminal_symbol.empty()) {
								set_is_empty = true;
								break;
							}
							set_for_nonterminal = sumOfSets(set_for_nonterminal, terminal_symbol, size); //iterative union of sets
						}
						if (!set_is_empty) {
							for (auto const& x : set_for_nonterminal) {
								first_k2[non_terminals[i]][x.first]=x.second;
							}
							
						}  //if both of the sets weren't empty we unite set on the previous step and on the current
					}
				}

			}
			if (first_k2 == temp_configuration) configuration = true; // we check whether the algorithm has stabilized
			else temp_configuration = first_k2;
		}
		first_k2.erase("eps");
	}

	//searching for epsilon non-terminals
	void epsilon_non_term() {
		vector<string> rule;
		set<string> state=epsilon; //to know when the algorithm stabilizes
		bool only_epsilons,is_stabilized=false;
		while (!is_stabilized) { //while not stabilized
			for (int i = 0; i < non_terminals.size(); ++i) { //for each non-terminal
				int k = rules[non_terminals[i]].size(); 
				for (int j = 0; j < rules[non_terminals[i]].size(); ++j) { //for each rule
					rule = rules[non_terminals[i]][j];
					only_epsilons = true;
					if (rule[0] == "eps") epsilon.insert(non_terminals[i]); //if we have non-terminal->eps it's an epsilon non-terminal by definition
					else for (int k = 0; k < rule.size(); ++k) {
						if (epsilon.count(rule[k]) == 0) {
							only_epsilons = false;
							break;
						}
					}
					if (only_epsilons == true) epsilon.insert(non_terminals[i]); //if the rule consists of only epsilons we can derive eps from the current non-terminal
				}
			}
			if (state == epsilon) is_stabilized = true; // we check whether the algorithm has stabilized
			else state = epsilon;
		}
	}

	//inserting reachable non-terminals as the first symbol of the selected non-terminal
	void insertNonTerminals(int i, set<string>& state1) { // set of reachable non-terminals as the first symbol of the selected non-terminal in leftRecursive_get() funct 
		int k = rules[non_terminals[i]].size(), pos;
		vector<string> rule;
		for (int j = 0; j < rules[non_terminals[i]].size(); ++j) { //for each rule
			rule = rules[non_terminals[i]][j];
			if (isNonTerminal(rule[0])) state1.insert(rule[0]); //if the first character in the rule is non-terminal then we add it to our set
			pos = 0;
			while ((pos < rule.size()) && (epsilon.count(rule[pos]) > 0) && isNonTerminal(rule[pos])) { //while non-terminals can derive eps adding them to our set
				state1.insert(rule[pos]);
				++pos;
			}
			if (pos < rule.size() && isNonTerminal(rule[pos])) state1.insert(rule[pos]); //the first non-terminal after a sequence of epsilon non-terminals is also reachable
		}
	}

	//getting left recursive non-terminals
	void leftRecursive_get() {
		set<string> state1,state2;
		vector<string> rule;
		int is_stabilized = 0, pos=0;
		for (int p = 0; p < non_terminals.size(); ++p) { //for each rule
			state1.clear();
			insertNonTerminals(p, state1);
			if (!state1.empty()) { //state1 is empty if there's no left recursion for sure
				is_stabilized = 0;
				while (!is_stabilized) { //while the algorithm isn't stable
					for (string x : state1) {
						for (int i = 0; i < non_terminals.size(); ++i) if (non_terminals[i] == x) {
							insertNonTerminals(i, state1); break; //for each non-terminal in our set we insert reachable non-terminlals from the first place in the rule
						}
						
					}
					if (state1 == state2) is_stabilized = true; // we check whether the algorithm has stabilized
					else state2 = state1;
				}
				if (state1.count(non_terminals[p]) > 0) leftRecursive.insert(non_terminals[p]); //if we can reach our selected non-terminal then it's left recursive
			}
		}
	}

	//pretty format functions
	void outWithSpaces(pair<string, vector<int>> x, bool coma) {
		int position = 0;
		int pos_mas = 0;
		while (pos_mas < x.second.size() - 1) {
			cout << x.first.substr(position, x.second[pos_mas]) << " ";
			position += x.second[pos_mas];
			++pos_mas;
		}
		if (coma) cout << x.first.substr(position, x.second[pos_mas]) << ", "; //do we need coma at the end
		else cout << x.first.substr(position, x.second[pos_mas]);
	}

	void outWithSpacesToFile(pair<string, vector<int>> x, bool coma, ofstream& file) {
		int position = 0;
		int pos_mas = 0;
		while (pos_mas < x.second.size() - 1) {
			file << x.first.substr(position, x.second[pos_mas]) << " ";
			position += x.second[pos_mas];
			++pos_mas;
		}
		if (coma) file << x.first.substr(position, x.second[pos_mas]) << ", "; //do we need coma at the end
		else file << x.first.substr(position, x.second[pos_mas]);
	}

	//first_k sets output
	void first_k_out(bool spaces) { //do we need spaces between terminals
		int size, position = 0, pos_mas = 0;;
		for (int i = 0; i < non_terminals.size(); ++i) {
			cout << non_terminals[i] << ": {";
			size = 0;
			for (auto const& x : first_k2[non_terminals[i]]) {
				
				if (size < first_k2[non_terminals[i]].size() - 1) {
					if (spaces)  outWithSpaces(x,true);
					else cout << x.first << ", ";
				}
				else if (spaces) outWithSpaces(x,false); else cout << x.first;
				++size;
			}

			cout << "}" << endl;
		}
	}

	//epsilon non-terminals set output to the file
	void epsilon_out() {
		cout << "epsilon non-terminals: {";
		int size = 0;
		for (string x : epsilon) {
			if (size < epsilon.size()-1) cout << x << ", ";
			else cout << x;
			++size;
		}
		cout << "}";
	}

	//first_k sets output to the file
	void first_k_out_file(string filename, bool spaces) { //do we need spaces between terminals
		ofstream file(filename);
		int size, position = 0, pos_mas = 0;
		for (int i = 0; i < non_terminals.size(); ++i) {
			file << non_terminals[i] << ": {";
			size = 0;
			for (auto const& x : first_k2[non_terminals[i]]) {

				if (size < first_k2[non_terminals[i]].size() - 1) {
					if (spaces)  outWithSpacesToFile(x, true, file);
					else file << x.first << ", ";
				}
				else if (spaces) outWithSpacesToFile(x, false, file); else file << x.first;
				++size;
			}

			file << "}" << endl;
		}
		file.close();
	}

	//epsilon non-terminals set output to the file
	void epsilon_out_file(string filename) {
		ofstream file(filename, ios::app);
		file << "epsilon non-terminals: {";
		int size = 0;
		for (string x : epsilon) {
			if (size < epsilon.size() - 1) file << x << ", ";
			else file << x;
			++size;
		}
		file << "}";
		file.close();
	}

	//left recursive set output
	void left_recursive_out() {
		cout << "left-recursive non-terminals: {";
		int size = 0;
		for (string x : leftRecursive) {
			if (size < leftRecursive.size() - 1) cout << x << ", ";
			else cout << x;
			++size;
		}
		cout << "}";
	}

	//left recursive set output to the file
	void left_recursive_out_file(string filename) {
		ofstream file(filename, ios::app);
		file << "left-recursive non-terminals: {";
		int size = 0;
		for (string x : leftRecursive) {
			if (size < leftRecursive.size() - 1) file << x << ", ";
			else file << x;
			++size;
		}
		file << "}";
		file.close();
	}
};

//output k to the file
void out_k(string filename, int k) {
	ofstream file(filename, ios::app);
	file << endl << "k = " << k <<endl;
	file.close();
}

int main() {
	int k = 3;
	Grammar grammar;
	grammar.read("Grammar.txt");
	grammar.build_first_k(k);
	grammar.epsilon_non_term();
	grammar.first_k_out_file("Output.txt",false);
	grammar.epsilon_out_file("Output.txt");
	out_k("Output.txt", k);
	grammar.leftRecursive_get();
	grammar.left_recursive_out_file("Output.txt");
}