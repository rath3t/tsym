
#include <sstream>
#include "abc.h"
#include "sum.h"
#include "baseptrlist.h"
#include "tsymtests.h"

using namespace tsym;

TEST_GROUP(BasePtrList) {};

TEST(BasePtrList, emptyListItemAccess)
{
    const BasePtrList empty;

    disableLog();
    CHECK(empty.front()->isUndefined());
    CHECK(empty.back()->isUndefined());
    enableLog();
}

TEST(BasePtrList, nonConstReverseIteration)
{
    const BasePtr plainList[] { a, b, c, d };
    BasePtrList list { d, c, b, a };
    int i = 0;

    for (auto it = list.rbegin(); it != list.rend(); ++it, ++i) {
        CHECK_EQUAL(plainList[i], *it);
    }
}

TEST(BasePtrList, popFrontAndBack)
{
    BasePtrList list { eight, b, c, two };

    CHECK_EQUAL(eight, list.pop_front());
    CHECK_EQUAL(3, list.size());
    CHECK_EQUAL(two, list.pop_back());
    CHECK_EQUAL(2, list.size());
    CHECK_EQUAL(c, list.pop_back());
    CHECK_EQUAL(1, list.size());
}

TEST(BasePtrList, restOfEmptyList)
{
    const BasePtrList empty;
    BasePtrList result;

    disableLog();
    result = empty.rest();
    enableLog();

    CHECK(result.empty());
}

TEST(BasePtrList, streamOperatorAndPushFront)
{
    const std::string expected("[ a   3 + b + c   2 ]");
    std::stringstream stream;
    BasePtrList list;

    list.push_front(two);
    list.push_front(Sum::create(b, c, three));
    list.push_front(a);

    stream << list;

    CHECK_EQUAL(expected, stream.str());
}

TEST(BasePtrList, constructWithInitializerList)
{
    BasePtrList list = { a, b, c, two, three, four };

    CHECK_EQUAL(6, list.size());
    CHECK_EQUAL(a, list.pop_front());
    CHECK_EQUAL(b, list.pop_front());
    CHECK_EQUAL(four, list.back());
}

TEST(BasePtrList, assignWithInitializerList)
{
    BasePtrList list(a, b);

    list = { c, d, two };

    CHECK_EQUAL(3, list.size());
    CHECK_EQUAL(c, list.pop_front());
    CHECK_EQUAL(d, list.pop_front());
    CHECK_EQUAL(two, list.front());
}
