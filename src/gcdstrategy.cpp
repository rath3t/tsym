
#include <cassert>
#include <cmath>
#include "gcdstrategy.h"
#include "numeric.h"
#include "logging.h"
#include "power.h"
#include "sum.h"
#include "polyinfo.h"
#include "product.h"

tsym::GcdStrategy::~GcdStrategy() {}

tsym::BasePtr tsym::GcdStrategy::compute(const BasePtr& u, const BasePtr& v) const
{
    PolyInfo polyInfo(u, v);

    if (polyInfo.isInputValid())
        return compute(u, v, polyInfo.listOfSymbols());

    TSYM_ERROR("Invalid gcd input: ", u, " and ", v);

    return Undefined::create();
}

tsym::BasePtr tsym::GcdStrategy::compute(const BasePtr& u, const BasePtr& v,
        const BasePtrList& L) const
{
    const BasePtr uExp(u->expand());
    const BasePtr vExp(v->expand());
    BasePtr result(Numeric::one());

    assert((L.empty() && uExp->isNumeric() && vExp->isNumeric()) || !L.empty());

    if (uExp->isOne() || vExp->isOne())
        ;
    else if (uExp->isZero() && vExp->isZero())
        return Undefined::create();
    else if (uExp->isZero())
        result = vExp;
    else if (vExp->isZero())
        result = uExp;
    else if (uExp->isEqual(vExp))
        result = uExp;
    else if (uExp->isNumeric() && vExp->isNumeric())
        result = computeNumerics(uExp, vExp);
    else if (!haveCommonSymbol(u, v, L))
        result = integerContent(u, v);
    else
        result = gcdViaAlgo(uExp, vExp, L);

    return normalize(result, L);
}

tsym::BasePtr tsym::GcdStrategy::computeNumerics(const BasePtr& u, const BasePtr& v) const
{
    const Number numU(u->numericEval());
    const Number numV(v->numericEval());
    Int intGcd(1);

    assert(numU.isRational() && numV.isRational());

    if (numU.isInt() && numV.isInt())
        intGcd = integerGcd(numU.numerator(), numV.numerator());

    return Numeric::create(intGcd);
}

tsym::Int tsym::GcdStrategy::integerGcd(const Int& u, const Int& v) const
{
    if (v == 0)
        return u;

    return integerGcd(v, u % v);
}

bool tsym::GcdStrategy::haveCommonSymbol(const BasePtr& u, const BasePtr& v,
        const BasePtrList& L) const
{
    for (const auto& item : L)
        if (u->has(item) && v->has(item))
            return true;

    return false;
}

tsym::BasePtr tsym::GcdStrategy::gcdViaAlgo(const BasePtr& u, const BasePtr& v,
        const BasePtrList& L) const
{
    const BasePtr intContent(integerContent(u, v));
    const BasePtr factor(Power::oneOver(intContent));
    const BasePtr uReduced(Product::create(factor, u));
    const BasePtr vReduced(Product::create(factor, v));
    BasePtr gcd;

    gcd = gcdAlgo(uReduced, vReduced, L);

    return Product::create(intContent, gcd);
}

tsym::BasePtr tsym::GcdStrategy::integerContent(const BasePtr& u, const BasePtr& v) const
{
    const Number uIntContent(integerContent(u));
    const Number vIntContent(integerContent(v));
    Int intGcd;

    if (!uIntContent.isInt() || !vIntContent.isInt())
        return Numeric::one();

    intGcd = integerGcd(uIntContent.numerator(), vIntContent.numerator());

    return Numeric::create(intGcd);
}

tsym::Number tsym::GcdStrategy::integerContent(const BasePtr& poly) const
{
    Number result;

    if (poly->isSum())
        result = integerContentOfSum(poly->operands());
    else
        result = poly->numericTerm()->numericEval().abs();

    return result.isInt() ? result : 1;
}

tsym::Number tsym::GcdStrategy::integerContentOfSum(const BasePtrList& summands) const
{
    Number intContent;
    Int result(0);

    for (const auto& summand : summands) {
        intContent = integerContent(summand);

        assert(intContent.isInt());

        result = integerGcd(result, intContent.numerator());
    }

    return Number(result);
}

tsym::BasePtr tsym::GcdStrategy::normalize(const BasePtr& result, const BasePtrList& L) const
{
    BasePtrList symbolListCopy(L);
    Number factor;

    factor = normalizationFactor(result, symbolListCopy);

    return Product::create(Numeric::create(factor), result);
}

tsym::Number tsym::GcdStrategy::normalizationFactor(const BasePtr& arg, BasePtrList& L) const
{
    BasePtr lCoeff;
    Number fac;

    if (L.empty())
        lCoeff = arg;
    else
        lCoeff = arg->leadingCoeff(L.pop_front());

    if (!lCoeff->isNumeric())
        return normalizationFactor(lCoeff, L);

    fac = lCoeff->numericEval();

    if (!fac.isRational())
        TSYM_ERROR(arg, " has a non-rational leading coefficient!");
    else if (fac.isZero())
        TSYM_ERROR("Gcd result has a zero leading coefficient: ", arg);

    return fac.sign();
}
