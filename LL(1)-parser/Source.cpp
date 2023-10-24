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
			   => I=0 => Strings={"+", "B", "A"}
			  |	
	Then S="A" 
			  |
			   => I=1 => Strings={"eps"}
	*/
	map<string, set<string>> first_k; //map for the first_k sets
	/*
	map<string, set<string>> where string is the name of non-terminal, set<string> is its first_k set
	*/
public:
	void read() {
		//reading axiom
		string axiom;
		string delimiter_axiom = "%start ";
		ifstream file("Grammar.txt"); //just hardcoded(
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
	}
	bool isNonTerminal(string s) {
		for (int i = 0; i < non_terminals.size(); ++i) if (s == non_terminals[i]) return true;
		return false;
	}
	bool isTerminal(string s) {
		for (int i = 0; i < terminals.size(); ++i) if (s == terminals[i]) return true;
		return false;
	}

	//sum binary operation on the words of a language, k is the length
	set<string> sumOfSets(set<string> set1, set<string> set2, int k) {
		set<string> res;
		for (string x : set1) {
			for (string y : set2) {
				if (x == "eps") res.insert(y.substr(0, k));
				else if (y == "eps") res.insert(x.substr(0, k));
				else res.insert((x + y).substr(0, k));
			}
		}
		return res;
	}
	// first_k building, size is the length
	void build_first_k(int size) {
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
						if (!set_is_empty) first_k[non_terminals[i]].insert(set_for_nonterminal.begin(), set_for_nonterminal.end());  //if both of the sets weren't empty we unite set on the previous step and on the current
					}
				}
				
			}
			if (first_k == temp_configuration) configuration = true; // we check whether the algorithm has stabilized
			else temp_configuration = first_k;
		}
	}
};
int main() {
	Grammar grammar;
	grammar.read();
	grammar.build_first_k(1);
}