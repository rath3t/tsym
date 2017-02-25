
#include <cstddef>
#include "base.h"
#include "baseptr.h"
#include "symbolregistry.h"
#include "symbolmap.h"
#include "undefined.h"
#include "printer.h"
#include "logging.h"

/* For use in the default constructor only. */
namespace tsym {
    namespace {
        const BasePtr& undefinedBaseForNoArgCtor()
        {
            static BasePtr undefined(Undefined::create());

            return undefined;
        }
    }
}

tsym::BasePtr::BasePtr() :
    bp(undefinedBaseForNoArgCtor().bp)
{
    /* This sets the reference count to 2, which is necessary to not destroy the static object
     * above, when this particular undefined BasePtr goes out of scope somewhere. */
    ++bp->refCount;

    #ifdef TSYM_DEBUG_STRINGS
    prettyStr = "Undefined";
    #endif
}

tsym::BasePtr::BasePtr(const Base *base) :
    bp(base)
{
    if (base == nullptr)
        TSYM_CRITICAL("Initiate Base class with null pointer!");

    ++bp->refCount;

    #ifdef TSYM_DEBUG_STRINGS
    prettyStr = (Printer(*this)).getStr();
    #endif

    if (bp->isSymbol())
        SymbolRegistry::add(base->name());
}

tsym::BasePtr::BasePtr(const BasePtr& other) :
    bp(other.bp)
{
    if (other.bp == nullptr)
        TSYM_CRITICAL("IniTIATE bASe class with null pointer!");

    #ifdef TSYM_DEBUG_STRINGS
    prettyStr = other.prettyStr;
    #endif

    ++bp->refCount;

    if (bp->isSymbol())
        SymbolRegistry::add(bp->name());
}

const tsym::BasePtr& tsym::BasePtr::operator = (const BasePtr& other)
{
    const Base* const old = bp;

    if (bp->isSymbol())
        SymbolRegistry::remove(bp->name());
    if (other->isSymbol())
        SymbolRegistry::add(other->name());

    bp = other.bp;

    #ifdef TSYM_DEBUG_STRINGS
    prettyStr = other.prettyStr;
    #endif

    ++bp->refCount;

    if (--old->refCount == 0)
        delete old;

    return *this;
}

tsym::BasePtr::~BasePtr()
{
    --bp->refCount;

    if (bp->isSymbol())
        SymbolRegistry::remove(bp->name());

    if (bp->refCount == 0)
        delete bp;
}

const tsym::Base *tsym::BasePtr::operator -> () const
{
    return bp;
}

const tsym::Base& tsym::BasePtr::operator* () const
{
    return *bp;
}

const tsym::Base *tsym::BasePtr::base() const
{
    return bp;
}

std::ostream& tsym::operator << (std::ostream& stream, const BasePtr& ptr)
{
    Printer printer(ptr);

    printer.print(stream);

    return stream;
}
