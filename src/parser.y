
%{
    #include <stdio.h>
    #include <string.h>
    #include "parseradapter.h"

    int yylex(void);
    int foundSyntaxError;
    void yyerror(char *message);
    void *result;
    extern char yytext[];
%}

%token LONG LONG_TEXT
%token DOUBLE
%token SYMBOL
%token DOUBLE_RANGE
%token SIN COS TAN ASIN ACOS ATAN ATAN2
%token SQRT
%token LOG
%token PI EULER

%union {
    char *stringVal;
    long longVal;
    double doubleVal;
    void *basePtr;
}

%type <stringVal> SYMBOL LONG_TEXT
%type <stringVal> DOUBLE_RANGE
%type <doubleVal> DOUBLE
%type <longVal> LONG
%type <basePtr> parsedExpr expr function textual

%left '+' '-'
%left '*' '/'
%nonassoc UNARY
%right '^'

%start parsedExpr

%%

parsedExpr: expr {
        result = $1;
    }
    ;

expr: function
    | textual
    | DOUBLE_RANGE {
        tsym_parserAdapter_logParsingError("Float out of representable range: ", $1);
        result = $$ = tsym_parserAdapter_createMaxDouble("Continue with maximum double: ");
    }
    | LONG {
        result = $$ = tsym_parserAdapter_createInteger($1);
    }
    | LONG_TEXT {
        result = $$ = tsym_parserAdapter_createLongInteger($1);
    }
    | DOUBLE {
        result = $$ = tsym_parserAdapter_createDouble($1);
    }
    | expr '+' expr {
        result = $$ = tsym_parserAdapter_createSum($1, $3);
        tsym_parserAdapter_deletePtrs($1, $3);
    }
    | expr '-' expr {
        result = $$ = tsym_parserAdapter_createDifference($1, $3);
        tsym_parserAdapter_deletePtrs($1, $3);
    }
    | expr '*' expr {
        result = $$ = tsym_parserAdapter_createProduct($1, $3);
        tsym_parserAdapter_deletePtrs($1, $3);
    }
    | expr '/' expr {
        result = $$ = tsym_parserAdapter_createQuotient($1, $3);
        tsym_parserAdapter_deletePtrs($1, $3);
    }
    | '+' expr %prec UNARY {
        result = $$ = $2;
    }
    | '-' expr %prec UNARY {
        result = $$ = tsym_parserAdapter_createMinus($2);
        tsym_parserAdapter_deletePtr($2);
    }
    | expr '^' expr {
        result = $$ = tsym_parserAdapter_createPower($1, $3);
        tsym_parserAdapter_deletePtrs($1, $3);
    }
    | '(' expr ')' {
        result = $$ = $2;
    }
    ;

textual: SYMBOL {
        result = $$ = tsym_parserAdapter_createSymbol($1);
    }
    | PI {
        result = $$ = tsym_parserAdapter_createPi();
    }
    | EULER {
        result = $$ = tsym_parserAdapter_createEuler();
    }
    ;

function: SIN '(' expr ')' {
        result = $$ = tsym_parserAdapter_createSine($3);
        tsym_parserAdapter_deletePtr($3);
    }
    | COS '(' expr ')' {
        result = $$ = tsym_parserAdapter_createCosine($3);
        tsym_parserAdapter_deletePtr($3);
    }
    | TAN '(' expr ')' {
        result = $$ = tsym_parserAdapter_createTangent($3);
        tsym_parserAdapter_deletePtr($3);
    }
    | ASIN '(' expr ')' {
        result = $$ = tsym_parserAdapter_createAsine($3);
        tsym_parserAdapter_deletePtr($3);
    }
    | ACOS '(' expr ')' {
        result = $$ = tsym_parserAdapter_createAcosine($3);
        tsym_parserAdapter_deletePtr($3);
    }
    | ATAN '(' expr ')' {
        result = $$ = tsym_parserAdapter_createAtangent($3);
        tsym_parserAdapter_deletePtr($3);
    }
    | ATAN2 '(' expr ',' expr ')' {
        result = $$ = tsym_parserAdapter_createAtangent2($3, $5);
        tsym_parserAdapter_deletePtr($3);
        tsym_parserAdapter_deletePtr($5);
    }
    | LOG '(' expr ')' {
        result = $$ = tsym_parserAdapter_createLogarithm($3);
        tsym_parserAdapter_deletePtr($3);
    }
    | SQRT '(' expr ')' {
        result = $$ = tsym_parserAdapter_createSquareRoot($3);
        tsym_parserAdapter_deletePtr($3);
    }
    | textual '(' expr ')' {
        foundSyntaxError = 1;
        tsym_parserAdapter_logParsingError("Unknown function, treated as variable!", "");
        tsym_parserAdapter_deletePtr($3);
    }
    ;

%%

void yyerror(char *s)
{
    foundSyntaxError = 1;

    if (strcmp(s, "syntax error") != 0)
        tsym_parserAdapter_logParsingError(s, yylval.stringVal);
    else
        tsym_parserAdapter_logParsingError("Yacc: found syntax error", "");
}
