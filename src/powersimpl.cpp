
#include <cmath>
#include <cassert>
#include "powersimpl.h"
#include "power.h"
#include "product.h"
#include "numeric.h"
#include "numpowersimpl.h"
#include "constant.h"
#include "logarithm.h"
#include "logging.h"

tsym::PowerSimpl::PowerSimpl() :
    one(Numeric::one())
{}

tsym::BasePtrList tsym::PowerSimpl::simplify(const BasePtr& base, const BasePtr& exp)
{
    if (doesInvolveComplexNumbers(base, exp))
        return BasePtrList(Undefined::create(), one);
    else if (base->isNumeric())
        return simplifyNumericBase(base, exp);
    else if (base->isPower())
        return simplifyPowerBase(base, exp);
    else if (base->isProduct())
        return simplifyProductBase(base, exp);
    else if (base->isUndefined())
        return BasePtrList(Undefined::create(), one);
    else if (base->isConstant())
        return simplifyConstantBase(base, exp);
    else
        /* No simplification applicable. */
        return BasePtrList(base, exp);
}

bool tsym::PowerSimpl::doesInvolveComplexNumbers(const BasePtr& base, const BasePtr& exp)
{
    if (base->isNegative() && exp->isNumericallyEvaluable())
        return !exp->numericEval().isInt();
    else
        return false;
}

tsym::BasePtrList tsym::PowerSimpl::simplifyNumericBase(const BasePtr& base, const BasePtr& exp)
{
    if (exp->isNumeric())
        return simplifyNumericPower(base, exp);
    else
        return BasePtrList(base, exp);
}

tsym::BasePtrList tsym::PowerSimpl::simplifyNumericPower(const BasePtr& base,
        const BasePtr& exp)
{
    const Number nBase(base->numericEval());
    const Number nExp(exp->numericEval());

    return simplifyNumericPower(nBase, nExp);
}

tsym::BasePtrList tsym::PowerSimpl::simplifyNumericPower(const Number& base, const Number& exp)
{
    NumPowerSimpl numericPow;
    BasePtr newBase;
    BasePtr newExp;
    BasePtr preFac;

    numericPow.setPower(base, exp);

    preFac = Numeric::create(numericPow.getPreFactor());
    newBase = Numeric::create(numericPow.getNewBase());
    newExp = Numeric::create(numericPow.getNewExp());

    if (preFac->isOne())
        return BasePtrList(newBase, newExp);
    else
        return BasePtrList(Product::create(preFac, Power::create(newBase, newExp)), one);
}

tsym::BasePtrList tsym::PowerSimpl::simplifyPowerBase(const BasePtr& powBase, const BasePtr& e2)
    /* Performs ((base)^e1)^e2 = (base)^(e1*e2) if it's possible. */
{
    const BasePtr e1(powBase->exp());
    const BasePtr base(powBase->base());
    BasePtr newBase;
    BasePtr newExp;

    if (doContractExpFirst(base, e1, e2))
        newBase = base;
    else if (doesChangeSign(base, e1))
        /* This enables ((-a - b)^2)^c = (a + b)^c with a, b > 0. */
        newBase = Product::minus(base);
    else if (doContractExpSecond(e1, e2))
        newBase = base;
    else
        /* No simplification possible. */
        return BasePtrList(powBase, e2);

    newExp = Product::create(e1, e2);

    if (newExp->isZero())
        /* For very small numeric exponents, this could be the case. */
        return BasePtrList(one, one);

    return simplify(newBase, newExp);
}

bool tsym::PowerSimpl::doContractExpFirst(const BasePtr& base, const BasePtr& e1, const BasePtr& e2)
{
    if (areTwoIntegerExp(e1, e2))
        return true;
    else if (areTwoFractionExpWithOddDenom(e1, e2))
        return true;
    else if (base->isPositive())
        return true;
    else
        return false;
}

bool tsym::PowerSimpl::areTwoIntegerExp(const BasePtr& exp1, const BasePtr& exp2)
{
    return isInteger(exp1) && isInteger(exp2);
}

bool tsym::PowerSimpl::isInteger(const BasePtr& arg)
{
    if (arg->isNumeric())
        return arg->numericEval().isInt();
    else
        return false;
}

bool tsym::PowerSimpl::areTwoFractionExpWithOddDenom(const BasePtr& exp1, const BasePtr& exp2)
{
    Number nExp1;
    Number nExp2;

    if (!exp1->isNumeric() || !exp2->isNumeric())
        return false;

    nExp1 = exp1->numericEval();
    nExp2 = exp2->numericEval();

    if (nExp1.isFrac() && nExp2.isFrac())
        return nExp1.denominator() % 2 != 0 && nExp2.denominator() % 2 != 0;
    return false;
}

bool tsym::PowerSimpl::doesChangeSign(const BasePtr& base, const BasePtr& exp1)
{
    Number nExp1;

    if (!base->isNegative())
        return false;
    else if (!exp1->isNumeric())
        return false;

    nExp1 = exp1->numericEval();

    return nExp1.isInt() && nExp1.numerator() % 2 == 0;
}

bool tsym::PowerSimpl::doContractExpSecond(const BasePtr& e1, const BasePtr& e2)
{
    const BasePtr newExp = Product::create(e1, e2);

    if (isOddInteger(e1))
        return newExp->isNumericallyEvaluable() && newExp->numericEval().abs() != 1;
    else if (isEvenInteger(e1))
        return isEvenInteger(newExp);
    else if (isFraction(e1) && (isInteger(e2) || isFraction(e2)))
        return isFraction(newExp);
    else if (e1->isNumericallyEvaluable() && e2->isNumericallyEvaluable())
        return !isEvenInteger(e1);
    else
        return false;
}

bool tsym::PowerSimpl::isOddInteger(const BasePtr& arg)
{
    return isInteger(arg) && arg->numericEval().numerator() % 2 != 0;
}

bool tsym::PowerSimpl::isEvenInteger(const BasePtr& arg)
{
    return isInteger(arg) && arg->numericEval().numerator() % 2 == 0;
}

bool tsym::PowerSimpl::isFraction(const BasePtr& arg)
{
    return arg->isNumeric() && arg->numericEval().isFrac();
}

tsym::BasePtrList tsym::PowerSimpl::simplifyProductBase(const BasePtr& base, const BasePtr& exp)
    /* Performs (a*b)^c = a^c*b^c if possible. */
{
    const bool doExpandAll = isInteger(exp);
    BasePtrList simplified;
    BasePtrList keep;

    for (const auto& factor : base->operands())
        if (doExpandAll || factor->isPositive())
            simplified.push_back(Power::create(factor, exp));
        else
            keep.push_back(factor);

    if (simplified.empty())
        return BasePtrList(base, exp);

    simplified.push_back(Power::create(Product::create(keep), exp));

    return BasePtrList(Product::create(simplified), one);
}

tsym::BasePtrList tsym::PowerSimpl::simplifyConstantBase(const BasePtr& base, const BasePtr& exp)
{
    if (isBaseEulerConstantAndExpLogarithm(base, exp))
        return BasePtrList(exp->operands().front(), one);
    else
        return BasePtrList(base, exp);
}

bool tsym::PowerSimpl::isBaseEulerConstantAndExpLogarithm(const BasePtr& base, const BasePtr& exp)
{
    if (base->isEqual(Constant::createE()) && exp->isFunction())
        return exp->isEqual(Logarithm::create(exp->operands().front()));
    else
        return false;
}
