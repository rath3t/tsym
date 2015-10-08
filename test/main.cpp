
#include "CppUTest/CommandLineTestRunner.h"
#include "printer.h"
#include "numeric.h"
#include "logarithm.h"
#include "constant.h"
#include "poly.h"
#include "name.h"
#include "abc.h"
#ifdef USE_TRLOG
#include "trlog/trlog.h"
#endif

void initConstructOnFirstUse()
    /* This initialization ensures that CppUTest doesn't spot a memory leak for the first test cases
     * that call the following functions employing local static variables. */
{
    tsym::BasePtr undefined;
    tsym::Number n(1);

    n /= n;
    n -= n;

    tsym::Name("dummy").tex();

    tsym::Numeric::zero();
    tsym::Numeric::one();
    tsym::Numeric::mOne();

    tsym::Logarithm::create(tsym::Constant::createE());

    tsym::poly::gcd(tsym::Numeric::one(), tsym::Numeric::one());
}

int main(int argc, char** argv)
{
    initConstructOnFirstUse();

    tsym::Printer::disableFractions();

#ifdef USE_TRLOG
    trlog::appendToCerr(trlog::WARNING, trlog::PRINT_LEVEL);
    trlog::appendToNewFile("doc/log/info.log", trlog::INFO);
    trlog::appendToNewFile("doc/log/debug.log", trlog::DEBUG);
#endif

    return CommandLineTestRunner::RunAllTests(argc, argv);
}
