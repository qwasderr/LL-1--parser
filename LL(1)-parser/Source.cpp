#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <algorithm>
using namespace std;

class T {
private:
    int number_of_rule;
    string non_terminal;
    map<string, vector<int>> L_number_of_rule;
public:
    T(int rule_number,  string non_terminal, map<string, vector<int>> L_number_of_rule) {
        this->number_of_rule = rule_number;
        this->non_terminal=std::move(non_terminal);
        this->L_number_of_rule = std::move(L_number_of_rule);
    }
    T() {
        this->number_of_rule = -1;
        this->non_terminal = "";
        this->L_number_of_rule["eps"] = {3};
    }

    [[nodiscard]] int getNumberOfRule() const {
        return number_of_rule;
    }

    void setNumberOfRule(int numberOfRule) {
        number_of_rule = numberOfRule;
    }

    [[nodiscard]] const string &getNonTerminal() const {
        return non_terminal;
    }

    void setNonTerminal(const string &nonTerminal) {
        non_terminal = nonTerminal;
    }

    [[nodiscard]] const map<string, vector<int>> &getLNumberOfRule() const {
        return L_number_of_rule;
    }

    void setLNumberOfRule(const map<string, vector<int>> &lNumberOfRule) {
        L_number_of_rule = lNumberOfRule;
    }

    // Comparison operator for '<'
    bool operator<(const T& other) const {
        if (this->number_of_rule != other.number_of_rule)
            return this->number_of_rule < other.number_of_rule;

        if (this->non_terminal != other.non_terminal)
            return this->non_terminal < other.non_terminal;

        return this->L_number_of_rule < other.L_number_of_rule;
    }

    // Comparison operator for '=='
    bool operator==(const T& other) const {
        return (this->number_of_rule == other.number_of_rule &&
                this->non_terminal == other.non_terminal &&
                this->L_number_of_rule == other.L_number_of_rule);
    }
};

class Cell {
public:
    Cell() {}

public:
    Cell(const vector<string> &values, int ruleNumber, const set<int> &positionOfTRules,
         const map<int, T> &positionOfTToReference) : values(values), rule_number(ruleNumber),
                                                      positions_of_t_rules(positionOfTRules),
                                                      position_of_T_to_reference(positionOfTToReference) {}

    [[nodiscard]] const vector<string> &getValues() const {
        return values;
    }

    void setValues(const vector<string> &val) {
        Cell::values = val;
    }

    [[nodiscard]] int getRuleNumber() const {
        return rule_number;
    }

    void setRuleNumber(int ruleNumber) {
        rule_number = ruleNumber;
    }

    [[nodiscard]] const set<int> &getPositionOfTRules() const {
        return positions_of_t_rules;
    }

    void setPositionOfTRules(const set<int> &positionOfTRules) {
        positions_of_t_rules = positionOfTRules;
    }

    [[nodiscard]] const map<int, T> &getPositionOfTToReference() const {
        return position_of_T_to_reference;
    }

    void setPositionOfTToReference(const map<int, T> &positionOfTToReference) {
        position_of_T_to_reference = positionOfTToReference;
    }

private:
    vector<string> values;
    int rule_number{};

    // position inside values
    set<int> positions_of_t_rules;
    map<int, T> position_of_T_to_reference;
};

class Grammar {
	vector <string> terminals;
	vector <string> non_terminals;
	string axiom;

    // map<left side of the rule(non-terminal), map<number of the rule for the non-terminal, vector<string> of terminals/non-terminals for this rule>>
	map<string, map<int, vector<string>>> rules;
    // rules with not local order, but general
	map<string, map<int, vector<string>>> rules_order;
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

    map <string, map<string, vector<int>>> follow_k;
    map <string, set<map<string, vector<int>>>> local_k;
    map<string, map <T, Cell>> table_of_control;
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


    void build_follow_k(int k) {

        for (auto &non_terminal: non_terminals) {
            follow_k[non_terminal] = {};
        }

        //We need to know for each rule on what min step it could be applied. See example of 4 iteration from book!
        map<string, int> min_step_for_non_terminals = calculate_min_requiered_steps();
        initialize_follow_k_with_axiom();

        bool changes_made;
        int iteration_number = 1;
        do {
            changes_made = false;

            //for comparing results of previous iteration with current.
            map<string, map<string, vector<int>>> follow_k_current;
            follow_k_current.insert(follow_k.begin(), follow_k.end());

            for (const auto &rule_entry: rules) {
                string rule_non_terminal = rule_entry.first;
                int min_step_for_start_non_terminal = min_step_for_non_terminals[rule_non_terminal];

                for (const auto &entry: rule_entry.second) {
                    int rule_number = entry.first;
                    vector<string> production = entry.second;

                    // For each production Non_terminal -> αBβ:
                    for (int i = 0; i < production.size(); i++) {
                        string maybe_non_terminal = production[i];
                        if (isNonTerminal(maybe_non_terminal)) {
                            // if on this step current Non_terminal could be produced, nothing to change. It's similar to "-" in the book. See example for C, D in the book.
                            if (iteration_number >= min_step_for_start_non_terminal) {

                                // retrieve continuation = β
                                vector<string> follow_of_non_terminal = continuation_of_production(production, i + 1);

                                // get first_k(β) = first_k(β_1) +_k first_k(β_2) +_k first_k(β_3) +_k ...
                                map<string, vector<int>> first_k_from_follow = get_first_k_for(follow_of_non_terminal, k);

                                map<string, vector<int>> follow_k_from_rule_non_terminal = follow_k_current[rule_non_terminal];
                                map<string, vector<int>> concat_result = sumOfSets(first_k_from_follow, follow_k_from_rule_non_terminal, k);


                                // Step 3: check should we stop to iterate or not
                                map<string, vector<int>> present_set_of_follow = follow_k_current[maybe_non_terminal];
                                map<string, vector<int>> mergedMap = mergeMaps(present_set_of_follow, concat_result);
                                if (mergedMap != present_set_of_follow) {
                                    follow_k[maybe_non_terminal] = mergedMap;
                                    changes_made = true;
                                }
                            }
                        }

                    }
                }
            }
            iteration_number += 1;
        } while (changes_made);
    }

    void initialize_follow_k_with_axiom() {
        map<string, vector<int>> &axiom_0_iteration = follow_k[axiom];
        axiom_0_iteration["eps"] = {3};
    }

    // get first_k(w) = first_k(w_1) +_k first_k(w_2) +_k first_k(w_3) +_k ...
    map<string, vector<int>> get_first_k_for(vector<string> w, int k) {
        map<string, vector<int>> first_k_for_w;

        // in the beginning just add first_k(w_1) to map
        first_k_for_w = mergeMaps(first_k_for_w, get_first_k_for_symbol(w[0]));
        // for each next w_j sum binary in with w_(j-1)
        for (int j = 1; j < w.size(); j++) {
            first_k_for_w = sumOfSets(first_k_for_w, get_first_k_for_symbol(w[j]), k);
        }
        
        return first_k_for_w;
    }

    map<string, vector<int>> get_first_k_for_symbol(const string &symbol) {
        map<string, vector<int>> symbol_set;
        if (isNonTerminal(symbol)) {
            symbol_set = first_k2[symbol];
        } else {
            symbol_set[symbol] = {static_cast<int>(symbol.length())};
        }
        return symbol_set;
    }

    static vector<string> continuation_of_production(vector<string> production, int start_from) {
        vector<string> continuation;
        for (int i = start_from; i < production.size(); i++) {
            continuation.push_back(production[i]);
        }
        if (continuation.empty()) continuation.emplace_back("eps");
        return continuation;
    }

    static void printMap(const map<string, vector<int>>& myMap) {
        for (const auto& entry : myMap) {
            cout << entry.first << ": ";

            for (int value : entry.second) {
                cout << value << " ";
            }

            cout << endl;
        }
    }

    static void print_3_step(map<string, vector<int>> mergedMap) {
        cout << "!!! NOT EQUAL" << endl;
        cout << "Updated map:" << endl;
        printMap(mergedMap);
    }

    static map<string, vector<int>> mergeMaps(const map<string, vector<int>>& map1, const map<string, vector<int>>& map2) {
        map<string, vector<int>> mergedMap = map1;

        for (const auto& entry : map2) {
            const string& key = entry.first;
            const vector<int>& values = entry.second;
            mergedMap[key] = values;
        }

        return mergedMap;
    }

    map<string, int> calculate_min_requiered_steps() {
        map<string, int> step_for_rule;
        step_for_rule[axiom] = 0;

        map<string, vector<string>> production_map;

        for (const auto& entry : rules) {
            const string& non_terminal = entry.first;
            const map<int, vector<string>>& rule_steps = entry.second;

            for (const auto& step_entry : rule_steps) {
                const vector<string>& productions = step_entry.second;

                for (const string& production : productions) {
                    if (isNonTerminal(production) && production != non_terminal) production_map[production].push_back(non_terminal);
                }
            }
        }

        for (const auto& entry : production_map) {
            const vector<string>& could_be_derived_from = entry.second;
            for (const string& non_terminal_from : could_be_derived_from) {
                if (non_terminal_from == axiom) {
                    step_for_rule[entry.first] = 2;
                }
            }
        }

        for (const auto& entry : production_map) {
            const string &non_terminal = entry.first;
            if(non_terminal != axiom) {
                int step = find_step_recursive(non_terminal, step_for_rule, production_map, {non_terminal});
                step_for_rule[non_terminal] = step;
            }
        }

        return step_for_rule;
    }

    int find_step_recursive(const string& non_terminal, map<string, int>& step_for_rule,
                            const map<string, vector<string>>& production_map, set<string> visited_non_terminals) {
        if (step_for_rule[non_terminal] > 1) {
            return step_for_rule[non_terminal];
        }

        const vector<string>& productions = production_map.at(non_terminal);
        int min_step = INT_MAX;

        for (const string& production : productions) {
            int step = 1;
            if(visited_non_terminals.find(production) != visited_non_terminals.end()) {
                visited_non_terminals.clear();
                return -1;
            }
            visited_non_terminals.insert(production);
            step = max(step, find_step_recursive(production, step_for_rule, production_map, visited_non_terminals) + 1);
            min_step = min(min_step, step);
        }

        return min_step;
    }

    /**
     * Algorithm:
     * σ_n(S, A_i) = σ_(n-1)(S, A_i) with L, where:
     *
     * A_j -> w_1 A_i w_2
     *
     * L = first_k(w_2) +_k L_p
     *
     * L_p ∈ Local_k(S, A_j)
     *
     * Hint: how to check yourself: - follow_k(A_i) = ∪ L_i, L_i ∈ Local_k(S, A_i)
     *
     * */
    void build_local_k(int k) {

        for (auto &non_terminal: non_terminals) {
            local_k[non_terminal] = {};
        }

        //We need to know for each rule on what min step it could be applied. 
        map<string, int> min_step_for_non_terminals = calculate_min_requiered_steps();
        initialize_local_k_with_axiom();

        bool changes_made;
        int iteration_number = 1;
        do {
            changes_made = false;

            //for comparing results of previous iteration with current.
            map <string, set<map<string, vector<int>>>> local_k_current;
            local_k_current.insert(local_k.begin(), local_k.end());

            for (const auto &rule_entry: rules) {
                string rule_non_terminal = rule_entry.first;
                int min_step_for_start_non_terminal = min_step_for_non_terminals[rule_non_terminal];

                for (const auto &entry: rule_entry.second) {
                    int rule_number = entry.first;
                    vector<string> production = entry.second;

                    // For each production Non_terminal -> αBβ:
                    for (int i = 0; i < production.size(); i++) {
                        string maybe_non_terminal = production[i];
                        if (isNonTerminal(maybe_non_terminal)) {
                            // if on this step current Non_terminal could be produced, nothing to change. It's similar to "-" in the book. See example for C, D in the book.
                            if (iteration_number >= min_step_for_start_non_terminal) {

                                // retrieve continuation = w
                                vector<string> follow_of_non_terminal = continuation_of_production(production, i + 1);

                                // get first_k(w) = first_k(w_1) +_k first_k(w_2) +_k first_k(w_3) +_k ...
                                map<string, vector<int>> first_k_from_w2 = get_first_k_for(follow_of_non_terminal, k);

                                // get set of L_p for A_j
                                set<map<string, vector<int>>> &local_k_current_for_rule_non_terminal = local_k_current[rule_non_terminal];

                                // set before adding:
                                set<map<string, vector<int>>> local_k_current_for_maybe_non_terminal = local_k_current[maybe_non_terminal];
                                // set where new sets will be added:
                                set<map<string, vector<int>>> &local_k_for_maybe_non_terminal = local_k[maybe_non_terminal];

                                for (const auto &set_of_possible_follows: local_k_current_for_rule_non_terminal) {
                                    const map<string, vector<int>> &new_set_of_possible_follows = sumOfSets(first_k_from_w2, set_of_possible_follows, k);
                                    local_k_for_maybe_non_terminal.insert(new_set_of_possible_follows);
                                }

                                // Step 3: check should we stop to iterate or not
                                if (local_k_for_maybe_non_terminal != local_k_current_for_maybe_non_terminal) {
                                    changes_made = true;
                                }
                            }
                        }
                    }
                }
            }
            iteration_number += 1;
        } while (changes_made);
    }

    void initialize_local_k_with_axiom() {
        set<map<string, vector<int>>> &axiom_0_iteration = local_k[axiom];
        map<string, vector<int>> map_with_eps;
        map_with_eps["eps"] = {3};
        axiom_0_iteration.insert(map_with_eps);
    }

    void build_table_of_control(int k) {
        // TODO: generateCombinations(terminals, k);
        // Generate  Σ_k
        vector<string> combinations = get_combination_of_non_terminal(terminals, k);
        initialize_table_of_control(combinations, k);

        int number_of_tables = calculate_number_of_tables();
        map<int, T> T_rules = getT_Rules(number_of_tables);
        int counter_of_t_rules = 1;

        // Iterate over each T_rule, starting from 0 to next ones.
        for (const auto &t_rule: T_rules) {

            string non_terminal = t_rule.second.getNonTerminal();
            map<string, vector<int>> L_number_of_rule = t_rule.second.getLNumberOfRule();
            map<int, vector<string>> rules_of_non_terminal = rules_order[non_terminal];

            // A -> w, where w = A_1 a_1 A_2 a_2 A_3 a_3 ...
            for (const auto &rule: rules_of_non_terminal) {
                // w
                vector<string> w = rule.second;

                // u = First_k(w) +_k L_previous
                map<string, vector<int>> u_as_map = sumOfSets(get_first_k_for(w, k), L_number_of_rule, k);
                vector<string> u = extract_keys_from_map(u_as_map);

                // position inside w -> T rule
                map<int, T> position_of_t_rules;
                set<int> set_of_positions;

                for (int i =0; i < w.size(); i++) {
                    string w_i = w[i];
                    if(isNonTerminal(w_i)) {

                        // retrieve continuation = w_(i+1_ w_(i+2) ...
                        vector<string> follow_of_non_terminal = continuation_of_production(w, i + 1);
                        map<string, vector<int>> first_k_for_continuation = get_first_k_for(follow_of_non_terminal, k);
                        // create new L for new T
                        const map<string, vector<int>> L_new = sumOfSets(first_k_for_continuation, L_number_of_rule, k);

                        // create T. If new rule requires index > max, it means we already check this (T4 = T2 in example, try to find correspond T)
                        T T_new = *new T(counter_of_t_rules, w_i, L_new);
                        if (counter_of_t_rules >= number_of_tables) {
                            for (const auto &item: T_rules) {
                                if(item.second.getNonTerminal() == w_i && item.second.getLNumberOfRule() == L_new) {
                                    T_new = item.second;
                                }
                            }
                        }
                        // update data for cell value
                        position_of_t_rules[i] = T_new;
                        set_of_positions.insert(i);
                        T_rules[counter_of_t_rules++] = T_new;
                    }
                }

                Cell new_cell = *new Cell(w, rule.first, set_of_positions, position_of_t_rules);
                // Update table of control for current rule. u_i ∈ u
                for (const auto &u_i: u) {
                    table_of_control[u_i][t_rule.second] = new_cell;
                }
            }

            // To stop traverse
            if(counter_of_t_rules > number_of_tables) break;
        }
    }

    vector<string> get_combination_of_non_terminal(vector<string> input, int k) {
        set<string> result(input.begin(), input.end());
        result.insert("eps");

        int iteration = 2;
        while (iteration <= k) {
            set<string> intermediate_result;
            for (const auto &w_i: result) {
                for (const auto &w_j: result) {
                    map<string, vector<int>> map_w_i;
                    map_w_i[w_i] = {static_cast<int>(w_i.length())};
                    map<string, vector<int>> map_w_j;
                    map_w_j[w_j] = {static_cast<int>(w_j.length())};
                    const map<string, vector<int>> concat = sumOfSets(map_w_i, map_w_j, iteration);
                    const vector<string> combinations = extract_keys_from_map(concat);
                    for (const auto &item: combinations) {
                        intermediate_result.insert(item);
                    }
                }
            }
            result.insert(intermediate_result.begin(), intermediate_result.end());
            iteration++;
        }

        return vector(result.begin(), result.end());
    }

    static string vectorToString(const vector<string>& vect) {
        string result;
        for (const auto &item: vect) {
            result += item + " ";
        }
        return result;
    }

    static vector<string> extract_keys_from_map(const map<string, vector<int>>& map) {
        vector<string> result;
        for (const auto &item: map) {
            result.push_back(item.first);
        }
        return result;
    }

    map<int, T> getT_Rules(int number_of_tables) {
        map<int, T> T_rules;
        map<string, vector<int>> empty_set;
        empty_set["eps"] = {3};
        T T_0 = *new T(0, axiom, empty_set);
        T_rules[0] = (T_0);
        for(int i = 1; i < number_of_tables; i++) {
            // fill out with empty values
            T T_i = *new T(i, "", empty_set);
            T_rules[i] = T_i;

        }
        return T_rules;
    }

    void build_order_rules() {
        int counter = 1;

        // first change for axiom
        map<int, vector<string>> axiom_rules = rules[axiom];
        map<int, vector<string>> axiom_with_correct_order;
        for (const auto &item: axiom_rules) {
            axiom_with_correct_order[counter++] = item.second;
        }
        rules_order[axiom] = axiom_with_correct_order;


        for (const auto &rule: rules) {
            string non_terminal_rule = rule.first;
            if(non_terminal_rule != axiom) {
                map<int, vector<string>> rules_non_order = rule.second;
                map<int, vector<string>> rules_correct_order;
                for (const auto &item: rules_non_order) {
                    rules_correct_order[counter++] = item.second;
                }

                rules_order[non_terminal_rule] = rules_correct_order;
            }
        }
    }

    void order_rules_out_file(const string& filename, bool spaces) {
        ofstream file(filename, ios::app);
        file << "Order rules for our grammar: " << endl;
        map<int, string> result;
        for (const auto &non_terminal_to_rules: rules_order) {
            string non_terminal = non_terminal_to_rules.first;
            map<int, vector<string>> ruless = non_terminal_to_rules.second;
            for (const auto &rule: ruless) {
                string pretty_string = non_terminal + " -> " + vectorToString(rule.second) + "\t\t" + "(" + to_string(rule.first) + ")";
                result[rule.first] = pretty_string;
            }
        }

        for (const auto &item: result) {
            file << item.second << endl;
        }

        file << endl;
        file.close();
    }

    void initialize_table_of_control(const vector<string>& combinations, int k) {
        for (auto& combination: combinations) {
            table_of_control[combination] = {};
        }
    }

    int calculate_number_of_tables() {
        int counter = 0;
        for (const auto &local_k_item: local_k) {
            counter += local_k_item.second.size();
        }
        return counter;
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
		int size;
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
        file << "First_k for no terminals:" << endl;
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
        file << endl;
		file.close();
	}

    //follow_k sets output to the file
    void table_of_control_out_file(const string& filename, bool spaces) { //do we need spaces between terminals
        ofstream file(filename, ios::app);

        file << "Table of control for LL(k) k > 1" << endl;

        for (const auto& entry : table_of_control) {
            const string combination = entry.first;
            const map<T, Cell> innerMap = entry.second;

            file << "Combination: " << combination << endl;

            for (const auto& innerEntry : innerMap) {
                const T& tObject = innerEntry.first;
                const Cell& cellObject = innerEntry.second;

                file << "T_" << tObject.getNumberOfRule() << ": ";
                for (int i = 0; i < cellObject.getValues().size(); i++) {
                    string w_i = cellObject.getValues()[i];
                    if (cellObject.getPositionOfTRules().count(i) > 0) {
                        map<int, T> position_to_T_rule = cellObject.getPositionOfTToReference();
                        file << "T_" << position_to_T_rule[i].getNumberOfRule() << " ";
                    } else {
                        file << w_i << " ";
                    }
                }
                file << " ; rule - " << cellObject.getRuleNumber();
                file << endl;
            }

            file << endl;
        }

        file.close();
    }

    //local_k sets output to the file
    void local_k_out_file(const string& filename, bool spaces) { //do we need spaces between terminals
        ofstream file(filename, ios::app);
        int size, position = 0, pos_mas = 0;
        file << "Local_k for no terminals:" << endl;
        for (int i = 0; i < non_terminals.size(); ++i) {
            file << non_terminals[i] << ": {";
            string set_to_out;
            for (const auto &set: local_k[non_terminals[i]]) {
                size = 0;
                set_to_out += " { ";
                for (const auto &x: set) {
                    if (size < set.size() - 1) {
                        if (spaces)  outWithSpacesToFile(x, true, file);
                        else set_to_out += x.first + ", ";
                    }
                    else if (spaces) outWithSpacesToFile(x, false, file); else set_to_out += x.first;
                    ++size;
                }
                set_to_out += " },";
            }

            if (!set_to_out.empty()) {
                set_to_out.pop_back();
            }

            file << set_to_out;

            file << " }" << endl;
        }
        file << endl;
        file.close();
    }

	//follow_k sets output to the file
	void follow_k_out_file(const string& filename, bool spaces) { //do we need spaces between terminals
		ofstream file(filename, ios::app);
		int size, position = 0, pos_mas = 0;
        file << "Follow_k for no terminals:" << endl;
		for (int i = 0; i < non_terminals.size(); ++i) {
			file << non_terminals[i] << ": {";
			size = 0;
			for (auto const& x : follow_k[non_terminals[i]]) {

				if (size < follow_k[non_terminals[i]].size() - 1) {
					if (spaces)  outWithSpacesToFile(x, true, file);
					else file << x.first << ", ";
				}
				else if (spaces) outWithSpacesToFile(x, false, file); else file << x.first;
				++size;
			}

			file << "}" << endl;
		}
        file << endl;
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
		file << "}" << endl;
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
	int k = 2;
    string grammar_file = "Grammar.txt";
    string output_file = "Output.txt";
	Grammar grammar;
	grammar.read(grammar_file);
    grammar.build_order_rules();
	grammar.build_first_k(k);
    grammar.first_k_out_file(output_file,false);
    grammar.build_follow_k(k);
    grammar.follow_k_out_file(output_file,false);
    grammar.build_local_k(k);
    grammar.local_k_out_file(output_file, false);
    grammar.order_rules_out_file(output_file, false);
    grammar.build_table_of_control(k);
    grammar.table_of_control_out_file(output_file, false);
    grammar.epsilon_non_term();
    grammar.epsilon_out_file(output_file);
    out_k(output_file, k);
	grammar.leftRecursive_get();
	grammar.left_recursive_out_file(output_file);
}