
#include "symbol.h"
#include "tsymtests.h"

using namespace tsym;

TEST_GROUP(Symbol) {};

TEST(Symbol, creationByString)
{
    const BasePtr ptr = Symbol::create("abc");
    const Name expected("abc");

    CHECK(ptr->isSymbol());
    CHECK_EQUAL(expected, ptr->name());
}

TEST(Symbol, numericEvaluation)
{
    const BasePtr ptr = Symbol::create("dummy");

    CHECK_FALSE(ptr->isNumericallyEvaluable());
}

TEST(Symbol, typeString)
{
    const std::string expected("Symbol");
    const BasePtr basePtr = Symbol::create("dummy");

    CHECK_EQUAL(expected, basePtr->typeStr());
}

TEST(Symbol, emptyName)
{
    const BasePtr undefined = Symbol::create("");

    CHECK(undefined->isUndefined());
}
