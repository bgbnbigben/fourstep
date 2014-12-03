#include "parser.h"
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include "shunting-yard.h"
#include <cctype>
#include <cstring>

enum class Variable {
    Free,
    Positive,
    Integer,
    Binary,
};

std::map<std::string, point_type> varMap;
std::vector<std::string> varOrder;
std::vector<std::string> objFunc;

bool isIntrinsic(char c) {
    switch (c) {
        case '+':
        case '-':
        case '*':
        case '/':
        case '.':
        case '^':
        case '(':
        case ')':
        case '<':
        case '>':
        case '=':
            return true;
    }
    return false;
}

bool isFunc(std::string s) {
    if (s.compare("exp") == 0) return true;
    if (s.compare("log") == 0) return true;
    return false;
}


void show_error(Status status) {
    std::string message;
    switch (status) {
        case ERROR_SYNTAX:
            message = "Syntax error";
            break;
        case ERROR_OPEN_PARENTHESIS:
            message = "Missing parenthesis";
            break;
        case ERROR_CLOSE_PARENTHESIS:
            message = "Extra parenthesis";
            break;
        case ERROR_UNRECOGNIZED:
            message = "Unknown character";
            break;
        case ERROR_NO_INPUT:
            message = "Empty expression";
            break;
        case ERROR_UNDEFINED_FUNCTION:
            message = "Unknown function";
            break;
        case ERROR_FUNCTION_ARGUMENTS:
            message = "Missing function arguments";
            break;
        case ERROR_UNDEFINED_CONSTANT:
            message = "Unknown constant";
            break;
        default:
            message = "Unknown error";
    }
    std::cerr << message << std::endl;
}

double eval(std::string s) {
    double result = 0.0;
    Status status = shunting_yard(s.c_str(), &result);
    if (status != OK) {
        show_error(status);
        std::cerr << s << std::endl;
        // TODO work out how to do this without killing the entire program
        std::exit(1);
    } else {
        return result;
    }
}

REAL_TYPE gamsFunc(const points_vector& x) {
    auto ret = 0.0;
    for (auto i = 0u; i < objFunc.size(); i++) {
        std::string replaced;
        for (auto j = 0u; j < objFunc[i].size(); j++) {
            if (isIntrinsic(objFunc[i][j])) {
                replaced.push_back(objFunc[i][j]);
            } else if (j + 3 < objFunc[i].size() && isFunc(objFunc[i].substr(j, 3))) {
                replaced += objFunc[i].substr(j, 3);
                j += 2;
            } else if (isdigit(objFunc[i][j])) {
                replaced.push_back(objFunc[i][j]);
            } else if (isspace(objFunc[i][j])) {
                continue;
            } else {
                /* Assume variable */
                auto start = j;
                while (!isIntrinsic(objFunc[i][j]) && !isspace(objFunc[i][j])) {
                    j++;
                }
                std::string var = objFunc[i].substr(start, j - start);
                auto good = true;
                match(x[std::find(varOrder.begin(), varOrder.end(), var) - varOrder.begin()], [&](Point<REAL_TYPE> p) {
                    good = std::isfinite(p());
                    replaced += std::to_string(p());
                }, [&](Point<DISCRETE_TYPE> p) {
                    good = std::isfinite(p());
                    replaced += std::to_string(p());
                });
                if (!good)
                    return std::numeric_limits<REAL_TYPE>::infinity();
                j--;
            }
        }
        if (i + 1 != objFunc.size()) {
            // This is a constraint
            size_t finder;
            size_t f2;
            if ((finder = replaced.find("==")) != std::string::npos) {
                // Penalise the living crap out of ==, since we want to *hurt*
                // binary failures.
                if (ret < std::numeric_limits<REAL_TYPE>::max() - std::pow(1 + eval(replaced.substr(0, finder)) - eval(replaced.substr(finder+2)), 5.)) {
                    ret += std::pow(1 + eval(replaced.substr(0, finder)) - eval(replaced.substr(finder+2)), 5.);
                } else {
                    return std::numeric_limits<REAL_TYPE>::infinity();
                }
            } else if ((finder = replaced.find("<")) != std::string::npos && (f2 = replaced.find("<", finder+1)) != std::string::npos) {
                bool eq = false;
                if (replaced[finder+1] == '=') {
                    /* <= */
                    eq = true;
                    auto left = eval(replaced.substr(0, finder));
                    auto right = eval(replaced.substr(finder+2, f2));
                    if (left > right) {
                        if (ret < std::numeric_limits<REAL_TYPE>::max() - std::pow(1 + left - right, 5.)) {
                            ret += std::pow(1 + left - right, 5.);
                        } else {
                            return std::numeric_limits<REAL_TYPE>::infinity();
                        }
                    }
                } else {
                    /* < */
                    auto left = eval(replaced.substr(0, finder));
                    auto right = eval(replaced.substr(finder+1, f2));
                    if (left >= right) {
                        if (ret < std::numeric_limits<REAL_TYPE>::max() - std::pow(1 + left - right, 5.) + 10.000000) {
                            ret += std::pow(1 + left - right, 5.) + 10.000000;
                        } else {
                            return std::numeric_limits<REAL_TYPE>::infinity();
                        }
                    }
                }

                if (replaced[f2+1] == '=') {
                    /* <= */
                    auto left = eval(replaced.substr(finder + 1 + eq, f2));
                    auto right = eval(replaced.substr(f2+2));
                    if (left > right) {
                        if (ret < std::numeric_limits<REAL_TYPE>::max() - std::pow(1 + left - right, 5.)) {
                            ret += std::pow(1 + left - right, 5.);
                        } else {
                            return std::numeric_limits<REAL_TYPE>::infinity();
                        }
                    }
                } else {
                    /* < */
                    auto left = eval(replaced.substr(finder + 1 + eq, f2));
                    auto right = eval(replaced.substr(f2+1));
                    if (left >= right) {
                        if (ret < std::numeric_limits<REAL_TYPE>::max() - std::pow(1 + left - right, 5.) + 10.000000) {
                            ret += std::pow(1 + left - right, 5.) + 10.000000;
                        } else {
                            return std::numeric_limits<REAL_TYPE>::infinity();
                        }
                    }
                }
            } else if ((finder = replaced.find("<")) != std::string::npos) {
                if (replaced.find("=") != std::string::npos) {
                    auto left = eval(replaced.substr(0, finder));
                    auto right = eval(replaced.substr(finder+2)); 
                    if (left > right) {
                        if (ret < std::numeric_limits<REAL_TYPE>::max() - std::pow(1 + left - right, 5.)) {
                            ret += std::pow(1 + left - right, 5.);
                        } else {
                            return std::numeric_limits<REAL_TYPE>::infinity();
                        }
                    }
                } else {
                    auto left = eval(replaced.substr(0, finder));
                    auto right = eval(replaced.substr(finder+1));
                    if (left >= right) {
                        if (ret < std::numeric_limits<REAL_TYPE>::max() - std::pow(1 + left - right, 5.) + 10.000000) {
                            ret += std::pow(1 + left - right, 5.) + 10.000000;
                        } else {
                            return std::numeric_limits<REAL_TYPE>::infinity();
                        }
                    }
                }
            } else if ((finder = replaced.find(">")) != std::string::npos) {
                if (replaced.find("=") != std::string::npos) {
                    auto left = eval(replaced.substr(0, finder));
                    auto right = eval(replaced.substr(finder+2));
                    if (left < right) {
                        if (ret < std::numeric_limits<REAL_TYPE>::max() - std::pow(1 + right - left, 5.)) {
                            ret += std::pow(1 + right - left, 5.);
                        } else {
                            return std::numeric_limits<REAL_TYPE>::infinity();
                        }
                    }
                } else {
                    auto left = eval(replaced.substr(0, finder));
                    auto right = eval(replaced.substr(finder+1));
                    if (left <= right) {
                        if (ret < std::numeric_limits<REAL_TYPE>::max() - std::pow(1 + right - left, 5.) + 10.000000) {
                            ret += std::pow(1 + right - left, 5.) + 10.000000;
                        } else {
                            return std::numeric_limits<REAL_TYPE>::infinity();
                        }
                    }
                }
            }
        } else {
            if (ret < std::numeric_limits<REAL_TYPE>::max() - eval(replaced)) {
                ret += eval(replaced);
            } else {
                return std::numeric_limits<REAL_TYPE>::infinity();
            }
        }
    }
    return ret;
}


bool starts_with(std::string prefix, std::string& s) {
    if (prefix.size() > s.size()) return false;
    return std::equal(prefix.begin(), prefix.end(), s.begin());
}

bool ends_with(std::string suffix, std::string& s) {
    if (suffix.size() > s.size()) return false;
    return std::equal(suffix.rbegin(), suffix.rend(), s.rbegin());
}

std::vector<std::string> split(const std::string &text, char sep) {
    std::vector<std::string> tokens;
    int start = 0, end = 0;
    while ((end = text.find(sep, start)) != std::string::npos) {
        tokens.push_back(text.substr(start, end - start));
        start = end + 1;
    }
    tokens.push_back(text.substr(start));
    return tokens;
}

void stripStart(std::string& in, std::string token) {
    while(in.compare(0,1,token)==0)
        in.erase(in.begin()); // remove leading whitespaces
}

points_vector parseGams(char* f) {
    std::ifstream file(f);
    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> tokens;
        Variable type;
        if (starts_with("//", line) || line.length() == 0) {
            continue;
        } else if (line.find("VARIABLES") != std::string::npos) {
            std::string data(line);
            while (!ends_with(";", data)) {
                std::getline(file, line);
                stripStart(line, " ");
                data.append(line);
            }
            data.pop_back();
            if (starts_with("BINARY_", data)) {
                tokens = split(data.substr(18), ',');
                type = Variable::Binary;
            } else if (starts_with("POSITIVE_", data)) {
                tokens = split(data.substr(20), ',');
                type = Variable::Positive;
            } else if (starts_with("INTEGER_", data)) {
                tokens = split(data.substr(19), ',');
                type = Variable::Integer;
            } else {
                tokens = split(data.substr(11), ',');
                type = Variable::Free;
            }

            for (auto token: tokens) {
                if (type == Variable::Binary || type == Variable::Integer) {
                    if (type == Variable::Binary) {
                        auto found = varMap.find(token);
                        if (found != varMap.end()) { 
                            assert(strcmp(found->second.type().name(), typeid(Point<DISCRETE_TYPE>(0, 0, 0)).name()) == 0);
                            auto old = boost::get<Point<DISCRETE_TYPE>>(found->second);
                            std::cout << " / done" << std::endl;
                            found->second = Point<DISCRETE_TYPE>(old(), 0, 1);
                        } else {
                            varMap.insert(std::make_pair(token, Point<DISCRETE_TYPE>(0, 0, 1)));
                            varOrder.push_back(token);
                        }
                    } else {
                        auto found = varMap.find(token);
                        if (found != varMap.end()) { 
                            // jfc
                            // this line is 400% bullshit.
                            assert(strcmp(found->second.type().name(), typeid(Point<DISCRETE_TYPE>(0, 0, 0)).name()) == 0);
                            auto old = boost::get<Point<DISCRETE_TYPE>>(found->second);
                            std::cout << " / done" << std::endl;
                            found->second = Point<DISCRETE_TYPE>(old(), old.left, old.right);
                        } else {
                            varMap.insert(std::make_pair(token, Point<DISCRETE_TYPE>(0, -100000, 100000)));
                            varOrder.push_back(token);
                        }
                    }
                } else if (type == Variable::Positive) {
                    auto found = varMap.find(token);
                    if (found != varMap.end()) { 
                        match(found->second, [&](Point<DISCRETE_TYPE> p) {
                            found->second = Point<DISCRETE_TYPE>(0, 0, p.right);
                        }, [&](Point<REAL_TYPE> p) {
                            found->second = Point<REAL_TYPE>(0, 0., p.right);
                        });
                    } else {
                        varMap.insert(std::make_pair(token, Point<REAL_TYPE>(0, 0, 100000)));
                        varOrder.push_back(token);
                    }
                } else {
                    varMap.insert(std::make_pair(token, Point<REAL_TYPE>(0, -100000, 100000)));
                    varOrder.push_back(token);
                }
            }
        } else if (starts_with("LOWER_BOUND", line)) {
            while (true) {
                std::getline(file, line);
                if (ends_with("}", line)) break;
                auto tokens = split(line, ':');
                assert(tokens.size() == 2);
                stripStart(tokens[1], " ");
                while(tokens[1].size()>0 && (tokens[1].compare(tokens[1].size()-1,1," ")==0 || tokens[1].compare(tokens[1].size()-1,1,";")==0))
                    tokens[1].erase(tokens[1].end()-1); // remove trailing whitespaces 
                auto f = varMap.find(tokens[0]);
                match(f->second, [&](Point<REAL_TYPE> p) {
                    //f->second = Point<REAL_TYPE>(std::max(p(), std::stod(tokens[1])), std::stod(tokens[1]), p.right);
                    f->second = Point<REAL_TYPE>(std::stod(tokens[1]), std::stod(tokens[1]), p.right);
                }, [&](Point<DISCRETE_TYPE> p) {
                    //f->second = Point<DISCRETE_TYPE>(std::max(p(), std::stoll(tokens[1])), std::stoll(tokens[1]), p.right);
                    f->second = Point<DISCRETE_TYPE>(std::stoll(tokens[1]), std::stoll(tokens[1]), p.right);
                });
            }
        } else if (starts_with("UPPER_BOUND", line)) {
            while (true) {
                std::getline(file, line);
                if (ends_with("}", line)) break;
                auto tokens = split(line, ':');
                assert(tokens.size() == 2);
                stripStart(tokens[1], " ");
                while(tokens[1].size()>0 && (tokens[1].compare(tokens[1].size()-1,1," ")==0 || tokens[1].compare(tokens[1].size()-1,1,";")==0))
                    tokens[1].erase(tokens[1].end()-1); // remove trailing whitespaces 
                auto f = varMap.find(tokens[0]);
                match(f->second, [&](Point<REAL_TYPE> p) {
                    f->second = Point<REAL_TYPE>(std::min(p(), std::stod(tokens[1])), p.left, std::stod(tokens[1]));
                }, [&](Point<DISCRETE_TYPE> p) {
                    f->second = Point<DISCRETE_TYPE>(std::min(p(), std::stoll(tokens[1])), p.left, std::stoll(tokens[1]));
                });
            }
        } else if (starts_with("EQUATION", line) || starts_with("OBJ", line)) {
            std::string data(line);
            while (!ends_with(";", data)) {
                std::getline(file, line);
                stripStart(line, " ");
                data.append(line);
            }
            data.pop_back();
            if (starts_with("EQUATION", data)) {
                auto rest = data.substr(9);
                stripStart(rest, " ");
                auto eqnTokens = split(rest, ',');
                std::vector<bool> eqns(eqnTokens.size(), false);
                for (auto todo = eqnTokens.size(); todo; todo--) {
                    std::string ll, eqnData;
                    while (ll.length() == 0 || !ends_with(";", eqnData)) {
                        std::getline(file, ll);
                        if (ll.length()) eqnData.append(ll);
                    }
                    eqnData.pop_back();
                    auto tokens = split(eqnData, ':');
                    stripStart(tokens[1], " ");
                    auto pos = 0; for (; pos < eqnTokens.size(); pos++) {
                        if (tokens[0].compare(eqnTokens[pos]) == 0) break;
                    }
                    eqns[pos] = true;
                    objFunc.push_back(tokens[1]);
                }
            } else {
                if (data.find("maximize") != std::string::npos) {
                    objFunc.push_back("-1 * (" + data.substr(14) + ")");
                } else {
                    objFunc.push_back(data.substr(14));
                }
            }
        } else {
            // Unimplemented token; ignore.
        }
    }
    /*for (auto it = varMap.begin(); it != varMap.end(); it++) {
        match (it->second, [&](Point<REAL_TYPE> p) {
            std::cout << it->first << "(real): " << p.left << "->" << p() << "->" << p.right << std::endl;
        }, [&](Point<DISCRETE_TYPE> p) {
            std::cout << it->first << "(disc): " << p.left << "->" << p() << "->" << p.right << std::endl;
        });
    }
    */

    points_vector x(varOrder.size(), Point<REAL_TYPE>(0, 0, 0));
    for (auto i = 0u; i < varOrder.size(); i++) {
        auto find = varMap.find(varOrder[i]);
        x[i] = find->second;
    }
    return x;

    /*
    std::cout << "This point gives " << gamsFunc(x) << std::endl;;
    */
}
