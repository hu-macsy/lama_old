
#include "tracing/CallTree.hpp"


#include "common/Exception.hpp"
#include "common/Walltime.hpp"

#include <fstream>
#include <stdint.h>
#include <unistd.h>
#include <vector>
#include <utility>

namespace tracing

{

using std::endl;

#define MAX_CT_COUNTERS 4

typedef uint64_t CounterArray[MAX_CT_COUNTERS];

static int traced_counters = 1;

static std::ofstream outfile;

static void printRegionEntry( int region_id, RegionEntry& region )
{
    outfile << "# begin info line" << endl;
    outfile << "fl " << region.getFileToken() << " " << region.getFileName() << endl;
    outfile << "fn " << region_id;
    outfile << " 0";  //  src_region_id ;
    outfile << " " <<  region.getRegionName();
    outfile << " 10";  // kind
    outfile << " " << region.getLine(); 
    outfile << " " << region.getLine(); 
    outfile << " ?"; 
    outfile << endl;
    outfile << "0";
    for ( int i = 0; i < traced_counters; ++i )
    {
        outfile << " 0";
    }
    outfile << endl;
    outfile << "# end info line" << endl;
}

static void printRegion( int region_id, bool callFlag )
{
    if ( callFlag )
    {
        outfile << "c";
    }

    outfile << "fl 0" << endl;;

    if ( callFlag )
    {
        outfile << "c";
    }

    outfile << "fn " << region_id << endl;
}

static void printCounters( int line, const CounterArray& counters )
{
    outfile <<  line;
   
    for( int i = 0; i < traced_counters; ++i )
    {
        outfile << " " << counters[i];
    }

    outfile << endl;
}

struct CTEntry
{
    int caller; // region id of the calling region
    int callee; // region id of the called region
    int scl;    // source code location for costs or for the call
    int calls;  // number of calls, as entries might be accumulated

    CounterArray costs;

    void writeEntry()
    {
        if ( callee == - 1 )
        {
             // cost line

             outfile << "# begin cost line" << endl;

             printRegion( caller, false );

             int scl = 0; // source code line not supported here

             printCounters( scl, costs );

             outfile << "# end cost line" << endl;
        }
        else
        {
             outfile << "# begin trace call" << endl;
             printRegion( caller, false );
             printRegion( callee, true );
             outfile << "calls " << calls << " 0" << endl; 
             int scl = 0; // source code line not supported here
             printCounters( scl, costs );
             outfile << "# end trace call" << endl;
        }
    }

    void addCall( int n, const CounterArray& vals )
    {
        for ( int i = 0; i < n; ++i )
        {
            costs[i] += vals[i];
        }

        calls++;  // add stands for exactly one additional call
    }

    void set( int caller, int callee, int scl, int n, const CounterArray& vals )
    {
        this->caller = caller;
        this->callee = callee;
        this->scl    = scl;

        for ( int i = 0; i < n; ++i )
        {
            costs[i] += vals[i];
        }

        calls = 1;
    }

    bool isSame( int other_caller, int other_callee, int other_scl )
    {
         return other_caller == caller && other_callee == callee && other_scl == scl;
    }
};

#define CALL_CACHE_SIZE 16

struct CallTreeInfo
{

private:

    CounterArray lastCounterVals;
    CounterArray totalCosts;

    CTEntry call_cache[CALL_CACHE_SIZE];
     
    int call_cache_pos;
    int call_cache_last;

    int cacheHit;
    int cacheMiss;

    int newPos()
    {
        if ( call_cache_pos < CALL_CACHE_SIZE ) 
        {
            return call_cache_pos++;
        }
 
        call_cache[call_cache_last].writeEntry();
 
        int pos = call_cache_last++;
   
        if ( call_cache_last == CALL_CACHE_SIZE )
        {
            call_cache_last = 0;
        }

        return pos;
    }

public:

    CallTreeInfo()
    {
        call_cache_pos = 0;
        call_cache_last = 0;
        cacheHit = 0;
        cacheMiss = 0;
    }

    void clear()
    {
        for ( int i = 0; i < call_cache_pos; ++i )
        {
            call_cache[i].writeEntry();
        }

        call_cache_pos  = 0;
        call_cache_last = 0;
    }

    int find( int caller, int callee, int scl )
    {
        for ( int i = 0; i < call_cache_pos; ++i )
        {
            if ( call_cache[i].isSame( caller, callee, scl ) )
            {
                cacheHit++;
                return i;
            }
        }

        cacheMiss++;
        return -1;
    }

    void add( int caller, int callee, int scl, const CounterArray& costs )
    {
        int pos = find( caller, callee, scl );
    
        if ( pos >= 0 )
        {
            call_cache[pos].addCall( traced_counters, costs );
        }

        else
        {
            pos = newPos();

            call_cache[pos].set( caller, callee, scl, traced_counters, costs );
        }
    }

    void getExclusiveCosts ( CounterArray& costs, const CounterArray& counterVals ) 
    { 
        for ( int i = 0; i < traced_counters; i++ )
        {
            costs[i] = counterVals[i] - lastCounterVals [i];
            totalCosts [i] += costs [i];
            lastCounterVals[i] = counterVals [i];
        }
    }

    void setCounters( const CounterArray& counterVals )
    {
        for ( int i = 0; i < traced_counters; i++ )
        {
            lastCounterVals [i] = counterVals [i];
        }
    }

    const CounterArray& getTotalCosts() const
    {
        return totalCosts;
    }

    const int getHits() const
    {
        return cacheHit;
    }

    const int getMisses() const
    {
        return cacheMiss;
    }
};

CallTreeInfo CT;

class Guard
{
public:
    Guard() {}

    ~Guard()
    {
        CT.clear();
        CallTree::finish();
    }
};

/** Class for call stack used for CallTree generation
 *
 *  An entry consists of the called region and the values of performance counters at begin of call
 */

class CT_CallStack
{

    typedef std::pair<int, CounterArray> StackEntryType;

    private:

    std::vector<StackEntryType> stack;

    public:

    void push( int region_id, const CounterArray& enterCounterVals )
    {
        StackEntryType entry;

        entry.first = region_id;

        stack.push_back( entry );
        
        CounterArray& enterVals = stack.back().second;

        for ( int i = 0; i < traced_counters; ++i )
        {
            enterVals[i] = enterCounterVals[i];
        }
    }

    void pop()
    {
        stack.pop_back();
    }

    /** Get the costs for the total call of the region  */

    void getCosts( CounterArray& costs, const CounterArray& exitCounterVals ) const
    {
        const CounterArray& enterVals = stack.back().second;

        for ( int i = 0; i < traced_counters; ++i )
        {
            costs[i] = exitCounterVals[i] - enterVals[i];
        }
    }

    /** Get the current region id, taken from last entry of call stack */

    int current() const
    {
        return stack.back().first;
    }

    /** Query if call stack is empty. */

    bool empty() const
    {
        return stack.empty();
    }
};

CT_CallStack myCallStack;

static void getCounters( uint64_t counter_vals[] )
{
    // currently we count only time ticks

    counter_vals[0] = common::Walltime::timestamp();
}

void CallTree::enter( const int region_id, RegionEntry& region )
{
    if ( outfile == NULL )
    {
        printf( "open calltree file\n" );

        outfile.open( "calltree.ct", std::ios::out );

        if ( outfile.fail() )
        {
            COMMON_THROWEXCEPTION( "Could not open calltree.ct" )
        }

        outfile << "pid " << getpid() << endl;
        outfile << "cmd program_name" << endl;
        outfile << "part 1" << endl;
        outfile << endl;

        outfile << "# rate " << common::Walltime::timerate() << endl;

        double rate = 1.0 / common::Walltime::timerate();

        outfile << "event WALL_TICKS wallticks" << endl;
        outfile << "events WALL_TICKS" << endl;
        outfile << "define WALL_TIME " << rate << " WALL_TICKS" << endl;
    }
  
    if ( region.firstAccess() )
    {
        // write Info about the region in the file
        printRegionEntry( region_id, region );
    }

    CounterArray counterVals;

    getCounters( counterVals );

    if ( myCallStack.empty() )
    {
        // This is the first time we enter a region, do initialization

        CT.setCounters( counterVals );
    }
    else
    {
        // Before we enter the called region add costs so far for caller region

        int caller_region = myCallStack.current();

        int scl = 4714;

        CounterArray costs;

        CT.getExclusiveCosts( costs, counterVals );

        CT.add( caller_region, -1, scl, costs );
    }

    myCallStack.push( region_id, counterVals );
}

void CallTree::leave( const int region_id, const RegionEntry& region )
{
    if ( outfile == NULL )
    {
        COMMON_THROWEXCEPTION( "calltree.ct not open" )
    }

    // Add a cost line for the current region 

    int scl = 0;

    CounterArray counterVals;
    CounterArray costs;

    getCounters( counterVals );

    CT.getExclusiveCosts( costs, counterVals );

    if ( myCallStack.current() != region_id )
    {
        COMMON_THROWEXCEPTION( "call stack mismatch: leave " << region_id << ", should be " << myCallStack.current() )
    }

    CT.add( region_id, -1, scl, costs );

    myCallStack.getCosts( costs, counterVals );

    myCallStack.pop();

    if ( !myCallStack.empty() )
    {
        CT.add( myCallStack.current(), region_id, scl, costs );
    }
}

void CallTree::finish()
{
    if ( outfile.is_open() )
    {
        const CounterArray& totalCosts = CT.getTotalCosts();

        outfile << "# closed by finish" << endl;
        outfile << "totals ";
        for ( int i = 0; i < traced_counters; ++i )
        {
            outfile << " " << totalCosts[i];
        }
        outfile << endl;
        outfile << "# cache has " << CT.getHits() << " hits and "
                << CT.getMisses() << " misses" << endl;
        outfile.close();
    }
}

static Guard myGuard;

}
