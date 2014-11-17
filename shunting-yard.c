// Copyright 2011 - 2014 Brian Marshall. All rights reserved.
//
// Use of this source code is governed by the BSD 2-Clause License that can be
// found in the LICENSE file.

#include "shunting-yard.h"
#include "stack.h"

#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef enum {
    TOKEN_NONE,
    TOKEN_OPEN_PARENTHESIS,
    TOKEN_CLOSE_PARENTHESIS,
    TOKEN_OPERATOR,
    TOKEN_NUMBER,
    TOKEN_IDENTIFIER
} TokenType;

typedef struct {
    TokenType type;
    char *value;
} Token;

typedef enum {
    OPERATOR_OTHER,
    OPERATOR_UNARY,
    OPERATOR_BINARY
} OperatorArity;

typedef enum {
    OPERATOR_NONE,
    OPERATOR_LEFT,
    OPERATOR_RIGHT
} OperatorAssociativity;

typedef struct {
    char symbol;
    int precedence;
    OperatorArity arity;
    OperatorAssociativity associativity;
} Operator;

static const Token NO_TOKEN = {TOKEN_NONE, NULL};

static const Operator OPERATORS[] = {
    {'!', 1, OPERATOR_UNARY,  OPERATOR_LEFT},
    {'^', 2, OPERATOR_BINARY, OPERATOR_RIGHT},
    {'+', 3, OPERATOR_UNARY,  OPERATOR_RIGHT},
    {'-', 3, OPERATOR_UNARY,  OPERATOR_RIGHT},
    {'*', 4, OPERATOR_BINARY, OPERATOR_LEFT},
    {'/', 4, OPERATOR_BINARY, OPERATOR_LEFT},
    {'%', 4, OPERATOR_BINARY, OPERATOR_LEFT},
    {'+', 5, OPERATOR_BINARY, OPERATOR_LEFT},
    {'-', 5, OPERATOR_BINARY, OPERATOR_LEFT},
    {'(', 6, OPERATOR_OTHER,  OPERATOR_NONE}
};

// Returns a list of tokens extracted from the expression. `token_count` will be
// set to the length of the list.
static Token *tokenize(const char *expression, int *token_count);

// Parses a tokenized expression.
static Status parse(const Token *tokens, int token_count, Stack **operands,
                    Stack **operators, Stack **functions);

// Pushes an operator to the stack after applying operators with a higher
// precedence.
static Status push_operator(const Operator *operator, Stack **operands,
                            Stack **operators);

// Pushes the multiplication operator to the stack.
static Status push_multiplication(Stack **operands, Stack **operators);

// Allocates memory for a double and pushes it to the stack.
static void push_double(double x, Stack **operands);

// Pops a double from the stack, frees its memory and returns its value.
static double pop_double(Stack **operands);

// Converts a string into a number and pushes it to the stack.
static Status push_number(const char *value, Stack **operands);

// Converts a constant identifier into its value and pushes it to the stack.
static Status push_constant(const char *value, Stack **operands);

// Applies an operator to the top one or two operands, depending on if the
// operator is unary or binary.
static Status apply_operator(const Operator *operator, Stack **operands);

// Applies a unary operator to the top operand.
static Status apply_unary_operator(const Operator *operator, Stack **operands);

// Applies a function to the top operand.
static Status apply_function(const char *function, Stack **operands);

// Returns the arity of an operator, using the previous token for context.
static OperatorArity get_arity(char symbol, const Token *previous);

// Returns a matching operator.
static const Operator *get_operator(char symbol, OperatorArity arity);

Status shunting_yard(const char *expression, double *result) {
    int token_count = 0;
    Token *tokens = tokenize(expression, &token_count);

    Stack *operands = NULL, *operators = NULL, *functions = NULL;
    Status status = parse(tokens, token_count, &operands, &operators,
                          &functions);
    if (operands)
        *result = round(pop_double(&operands) * 10e14) / 10e14;
    else if (status == OK)
        status = ERROR_NO_INPUT;

    for (int i = 0; i < token_count; i++)
        free(tokens[i].value);
    free(tokens);
    while (operands)
        pop_double(&operands);
    while (operators)
        stack_pop(&operators);
    while (functions)
        stack_pop(&functions);
    return status;
}

Token *tokenize(const char *expression, int *token_count) {
    Token *tokens = NULL;
    int count = 0;
    const char *c = expression;
    while (*c) {
        TokenType type = TOKEN_NONE;
        char *value = NULL;

        if (*c == '(')
            type = TOKEN_OPEN_PARENTHESIS;
        else if (*c == ')')
            type = TOKEN_CLOSE_PARENTHESIS;
        else if (strchr("!^*/%+-", *c)) {
            type = TOKEN_OPERATOR;
            value = calloc(2, sizeof *value);
            *value = *c;
        } else if (sscanf(c, "%m[0-9.]", &value))
            type = TOKEN_NUMBER;
        else if (sscanf(c, "%m[A-Za-z]", &value))
            type = TOKEN_IDENTIFIER;

        if (!isspace(*c)) {
            tokens = realloc(tokens, sizeof *tokens * ++count);
            tokens[count - 1].type = type;
            tokens[count - 1].value = value;
        }
        c += value ? strlen(value) : 1;
    }
    *token_count = count;
    return tokens;
}

Status parse(const Token *tokens, int token_count, Stack **operands,
             Stack **operators, Stack **functions) {
    Status status = OK;
    for (int i = 0; i < token_count; i++) {
        const Token *previous = i > 0 ? &tokens[i - 1] : &NO_TOKEN;
        const Token *next = i + 1 < token_count ? &tokens[i + 1] : &NO_TOKEN;
        switch (tokens[i].type) {
            case TOKEN_OPEN_PARENTHESIS:
                // Implicit multiplication: "(2)(2)".
                if (previous->type == TOKEN_CLOSE_PARENTHESIS)
                    status = push_multiplication(operands, operators);

                stack_push(operators, get_operator('(', OPERATOR_OTHER));
                break;
            case TOKEN_CLOSE_PARENTHESIS: {
                // Apply operators until the previous open parenthesis is found.
                bool found_parenthesis = false;
                while (*operators && status == OK && !found_parenthesis) {
                    const Operator *operator = stack_pop(operators);
                    if (operator->symbol == '(')
                        found_parenthesis = true;
                    else
                        status = apply_operator(operator, operands);
                }
                if (!found_parenthesis)
                    status = ERROR_CLOSE_PARENTHESIS;
                else if (*functions)
                    status = apply_function(stack_pop(functions), operands);
                break;
            }
            case TOKEN_OPERATOR:
                status = push_operator(
                    get_operator(*tokens[i].value,
                                 get_arity(*tokens[i].value, previous)),
                    operands, operators);
                break;
            case TOKEN_NUMBER:
                if (previous->type == TOKEN_CLOSE_PARENTHESIS ||
                        previous->type == TOKEN_NUMBER ||
                        previous->type == TOKEN_IDENTIFIER)
                    status = ERROR_SYNTAX;
                else {
                    status = push_number(tokens[i].value, operands);

                    // Implicit multiplication: "2(2)" or "2a".
                    if (next->type == TOKEN_OPEN_PARENTHESIS ||
                            next->type == TOKEN_IDENTIFIER)
                        status = push_multiplication(operands, operators);
                }
                break;
            case TOKEN_IDENTIFIER:
                // The identifier could be either a constant or function.
                status = push_constant(tokens[i].value, operands);
                if (status == ERROR_UNDEFINED_CONSTANT &&
                        next->type == TOKEN_OPEN_PARENTHESIS) {
                    stack_push(functions, tokens[i].value);
                    status = OK;
                } else if (next->type == TOKEN_OPEN_PARENTHESIS ||
                           next->type == TOKEN_IDENTIFIER) {
                    // Implicit multiplication: "a(2)" or "a b".
                    status = push_multiplication(operands, operators);
                }
                break;
            default:
                status = ERROR_UNRECOGNIZED;
        }
        if (status != OK)
            return status;
    }

    // Apply all remaining operators.
    while (*operators && status == OK) {
        const Operator *operator = stack_pop(operators);
        if (operator->symbol == '(')
            status = ERROR_OPEN_PARENTHESIS;
        else
            status = apply_operator(operator, operands);
    }
    return status;
}

Status push_operator(const Operator *operator, Stack **operands,
                     Stack **operators) {
    if (!operator)
        return ERROR_SYNTAX;

    Status status = OK;
    while (*operators && status == OK) {
        const Operator *stack_operator = stack_top(*operators);
        if (operator->arity == OPERATOR_UNARY ||
                operator->precedence < stack_operator->precedence ||
                (operator->associativity == OPERATOR_RIGHT &&
                 operator->precedence == stack_operator->precedence))
            break;

        status = apply_operator(stack_pop(operators), operands);
    }
    stack_push(operators, operator);
    return status;
}

Status push_multiplication(Stack **operands, Stack **operators) {
    return push_operator(get_operator('*', OPERATOR_BINARY), operands,
                         operators);
}

void push_double(double x, Stack **operands) {
    double *pointer = malloc(sizeof *pointer);
    *pointer = x;
    stack_push(operands, pointer);
}

double pop_double(Stack **operands) {
    const double *pointer = stack_pop(operands);
    double x = *pointer;
    free((void *)pointer);
    return x;
}

Status push_number(const char *value, Stack **operands) {
    char *end_pointer = NULL;
    double x = strtod(value, &end_pointer);

    // If not all of the value is converted, the rest is invalid.
    if (value + strlen(value) != end_pointer)
        return ERROR_SYNTAX;
    push_double(x, operands);
    return OK;
}

Status push_constant(const char *value, Stack **operands) {
    double x = 0.0;
    if (strcasecmp(value, "e") == 0)
        x = M_E;
    else if (strcasecmp(value, "pi") == 0)
        x = M_PI;
    else if (strcasecmp(value, "tau") == 0)
        x = M_PI * 2;
    else
        return ERROR_UNDEFINED_CONSTANT;
    push_double(x, operands);
    return OK;
}

Status apply_operator(const Operator *operator, Stack **operands) {
    if (!operator || !*operands)
        return ERROR_SYNTAX;
    if (operator->arity == OPERATOR_UNARY)
        return apply_unary_operator(operator, operands);

    double y = pop_double(operands);
    if (!*operands)
        return ERROR_SYNTAX;
    double x = pop_double(operands);
    Status status = OK;
    switch (operator->symbol) {
        case '^':
            x = pow(x, y);
            break;
        case '*':
            x = x * y;
            break;
        case '/':
            x = x / y;
            break;
        case '%':
            x = fmod(x, y);
            break;
        case '+':
            x = x + y;
            break;
        case '-':
            x = x - y;
            break;
        default:
            return ERROR_UNRECOGNIZED;
    }
    push_double(x, operands);
    return status;
}

Status apply_unary_operator(const Operator *operator, Stack **operands) {
    double x = pop_double(operands);
    switch (operator->symbol) {
        case '+':
            break;
        case '-':
            x = -x;
            break;
        case '!':
            x = tgamma(x + 1);
            break;
        default:
            return ERROR_UNRECOGNIZED;
    }
    push_double(x, operands);
    return OK;
}

Status apply_function(const char *function, Stack **operands) {
    if (!*operands)
        return ERROR_FUNCTION_ARGUMENTS;

    double x = pop_double(operands);
    if (strcasecmp(function, "abs") == 0)
        x = fabs(x);
    else if (strcasecmp(function, "sqrt") == 0)
        x = sqrt(x);
    else if (strcasecmp(function, "ln") == 0 ||
             strcasecmp(function, "log") == 0)
        x = log(x);
    else if (strcasecmp(function, "lb") == 0)
        x = log2(x);
    else if (strcasecmp(function, "lg") == 0)
        x = log10(x);
    else if (strcasecmp(function, "cos") == 0)
        x = cos(x);
    else if (strcasecmp(function, "sin") == 0)
        x = sin(x);
    else if (strcasecmp(function, "tan") == 0)
        x = tan(x);
    else if (strcasecmp(function, "exp") == 0)
        x = exp(x);
    else
        return ERROR_UNDEFINED_FUNCTION;
    push_double(x, operands);
    return OK;
}

OperatorArity get_arity(char symbol, const Token *previous) {
    if (symbol == '!' || previous->type == TOKEN_NONE ||
            previous->type == TOKEN_OPEN_PARENTHESIS ||
            (previous->type == TOKEN_OPERATOR && *previous->value != '!'))
        return OPERATOR_UNARY;
    return OPERATOR_BINARY;
}

const Operator *get_operator(char symbol, OperatorArity arity) {
    for (size_t i = 0; i < sizeof OPERATORS / sizeof OPERATORS[0]; i++) {
        if (OPERATORS[i].symbol == symbol && OPERATORS[i].arity == arity)
            return &OPERATORS[i];
    }
    return NULL;
}
