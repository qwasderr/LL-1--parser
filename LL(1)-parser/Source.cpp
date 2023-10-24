#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <string>
using namespace std;
class Grammar {
	vector <string> terminals;
	vector <string> non_termainals;
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
			this->non_termainals.push_back(non_terminal);
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
};
int main() {
	Grammar grammar;
	grammar.read();
}