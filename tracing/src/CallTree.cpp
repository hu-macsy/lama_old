
#include "tracing/CallTree.hpp"


#include "common/Exception.hpp"
#include "common/Walltime.hpp"

#include <cstdio>
#include <stdint.h>
#include <vector>
#include <utility>

namespace tracing

{

#define MAX_CT_COUNTERS 4

typedef uint64_t CounterArray[MAX_CT_COUNTERS];

static int traced_counters = 1;

static FILE* outfile = NULL;

static void printRegionEntry( int region_id, RegionEntry& region )
{
    fprintf( outfile, "# begin info line\n" );
    fprintf( outfile, "fl %d %s\n", region.getFileToken(), region.getFileName() );
    fprintf( outfile, "fn %d", region_id );

    fprintf( outfile, " %d", 0 ); // src_region_id
    fprintf( outfile, " %s", region.getRegionName() );
    fprintf( outfile, " %d", 10 ); // kind
    fprintf( outfile, " %d", region.getLine() ); 
    fprintf( outfile, " %d", region.getLine() ); 
    fprintf( outfile, " %s", "?" ); 
    fprintf( outfile, "\n" );
    fprintf( outfile, "0" );
    for ( int i = 0; i < traced_counters; ++i )
    {
        fprintf( outfile, " 0" );
    }
    fprintf( outfile, "\n" );
    fprintf( outfile, "# end info line\n" );
}

static void printRegion( int region_id, bool callFlag )
{
    if ( callFlag )
    {
        fprintf( outfile, "c" );
    }

    fprintf( outfile, "fl %d\n", 0 );

    if ( callFlag )
    {
        fprintf( outfile, "c" );
    }

    fprintf( outfile, "fn %d\n", region_id );
}

static void printCounters( int line, const CounterArray& counters )
{
    fprintf( outfile, "%d", line );
   
    for( int i = 0; i < traced_counters; ++i )
    {
        fprintf( outfile, " %lld", counters[i] );
    }

    fprintf( outfile, "\n" );
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

             fprintf( outfile, "# begin cost line\n" );

             printRegion( caller, false );

             int scl = 0; // source code line not supported here

             printCounters( scl, costs );

             fprintf( outfile, "# end cost line\n" );
        }
        else
        {
             fprintf( outfile, "# begin trace call\n" );
             printRegion( caller, false );
             printRegion( callee, true );
             fprintf( outfile, "calls %d 0\n", calls ); 
             int scl = 0; // source code line not supported here
             printCounters( scl, costs );
             fprintf( outfile, "# end trace call\n" );
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

        outfile = fopen( "calltree.ct", "w" );

        if ( outfile == NULL )
        {
            COMMON_THROWEXCEPTION( "Could not open calltree.ct" )
        }

        fprintf( outfile, "pid %d\n", getpid() );
        fprintf( outfile, "cmd program_name\n" );
        fprintf( outfile, "part %d\n", 1 );
        fprintf( outfile, "\n" );

        fprintf( outfile, "# rate %lld\n", common::Walltime::timerate() );

        double rate = 1.0 / common::Walltime::timerate();

        fprintf( outfile, "event WALL_TICKS wallticks\n" );
        fprintf( outfile, "events WALL_TICKS\n" );
        fprintf( outfile, "define WALL_TIME %f WALL_TICKS\n", rate );
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
    if ( outfile != NULL )
    {
        const CounterArray& totalCosts = CT.getTotalCosts();

        fprintf( outfile, "# closed by finish\n" );
        fprintf( outfile, "totals " );
        for ( int i = 0; i < traced_counters; ++i )
        {
            fprintf( outfile, " %lld", totalCosts[i] );
        }
        fprintf( outfile, "\n" );
        fprintf( outfile, "# cache has %d hits and %d misses\n", CT.getHits(), CT.getMisses() );
        fclose( outfile );
    }
 
    outfile = NULL;
}

static Guard myGuard;

}
