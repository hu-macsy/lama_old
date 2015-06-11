
#include "tracing/CallTree.hpp"


#include "common/Exception.hpp"

#include <cstdio>
#include <stdint.h>

namespace tracing

{

static int traced_counters = 1;

static FILE* outfile = NULL;

static void printRegionEntry( int region_id, RegionEntry& region )
{
    fprintf( outfile, "#begin info line\n" );
    fprintf( outfile, "fl %d %s\n", region.getFileToken(), region.getFileName() );
    fprintf( outfile, "fn %d", region_id );

    fprintf( outfile, " %d", 0 ); // src_region_id
    fprintf( outfile, " %s", region.getRegionName() );
    fprintf( outfile, " %d", 10 ); // kind
    fprintf( outfile, " %d", region.getLine() ); 
    fprintf( outfile, " %d", region.getLine() ); 
    fprintf( outfile, " %s", "?" ); 
    fprintf( outfile, "\n" );
    fprintf( outfile, "0 0 0\n" );
    fprintf( outfile, "#end info line\n" );
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

static void printCounters( int line, uint64_t counters[] )
{
    fprintf( outfile, "%d", line );
   
    for( int i = 0; i < traced_counters; ++i )
    {
        fprintf( outfile, " %lld", counters[i] );
    }

    fprintf( outfile, "\n" );
}

#define MAX_CT_COUNTERS 4

struct CTEntry
{
    int caller; // region id of the calling region
    int callee; // region id of the called region
    int scl;    // source code location for costs or for the call
    int calls;  // number of calls, as entries might be accumulated

    uint64_t counters_val[MAX_CT_COUNTERS];

    void writeEntry()
    {
        if ( callee == - 1 )
        {
             // cost line

             fprintf( outfile, "# begin cost line\n" );

             printRegion( caller, false );

             int scl = 0; // source code line not supported here

             printCounters( scl, counters_val );

             fprintf( outfile, "# end cost line\n" );
        }
        else
        {
             fprintf( outfile, "# begin trace call\n" );
             printRegion( caller, false );
             printRegion( callee, true );
             fprintf( outfile, "calls %d 0\n", calls );
             int scl = 0; // source code line not supported here
             printCounters( scl, counters_val );
             fprintf( outfile, "# end trace call\n" );
        }
    }

    void addCall( int n, uint64_t vals[] )
    {
        for ( int i = 0; i < n; ++i )
        {
            counters_val[i] += vals[i];
        }

        calls++;  // add stands for exactly one additional call
    }

    void set( int caller, int callee, int scl, int n, uint64_t vals[] )
    {
        this->caller = caller;
        this->callee = callee;
        this->scl = scl;

        for ( int i = 0; i < n; ++i )
        {
            counters_val[i] += vals[i];
        }

        calls = 0;
    }

    bool isSame( int other_caller, int other_callee )
    {
         return other_caller == caller && other_callee == callee;
    }
};

#define CALL_CACHE_SIZE 16

struct CallTreeInfo
{
    uint64_t last_counter_vals[ MAX_CT_COUNTERS ];
    uint64_t total_counter_vals[ MAX_CT_COUNTERS ];
  
    CTEntry call_cache[CALL_CACHE_SIZE];
     
    int call_cache_pos;
    int call_cache_last;

    int call_cache_hit;
    int call_cache_miss;

    CallTreeInfo()
    {
        call_cache_pos = 0;
        call_cache_last = 0;
        call_cache_hit = 0;
        call_cache_miss = 0;
    }

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

    void clear()
    {
        for ( int i = 0; i < call_cache_pos; ++i )
        {
            call_cache[i].writeEntry();
        }
    }

    int find( int caller, int callee )
    {
        for ( int i = 0; i < call_cache_pos; ++i )
        {
            if ( call_cache[i].isSame( caller, callee ) )
            {
                call_cache_hit++;
                return i;
            }
        }

        call_cache_miss++;
        return -1;
    }

    void add( int caller, int callee, int scl, uint64_t vals[] )
    {
        int pos = find( caller, callee );
    
        if ( pos >= 0 )
        {
            call_cache[pos].addCall( traced_counters, vals );
        }

        else
        {
            pos = newPos();

            call_cache[pos].set( caller, callee, scl, traced_counters, vals );
        }
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

std::vector<int> CallTree::theCallStack;

void CallTree::enter( const int region_id, RegionEntry& region, double time )
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

        fprintf( outfile, "event WALL_TICKS wallticks\n" );
    }
  
    if ( region.firstAccess() )
    {
        // write Info about the region in the file
        printRegionEntry( region_id, region );
    }

    if ( theCallStack.size() == 0 )
    {
        // This is the first time we enter a region, we could do some initialization
    }
    else
    {
        int call_region = theCallStack.back();

        uint64_t vals[2] = { 10, 10 };
        int scl = 4714;

        CT.add( call_region, region_id, scl, vals );
    }

    theCallStack.push_back( region_id );
}

void CallTree::leave( const int region_id, const RegionEntry& region, double time )
{
    if ( outfile == NULL )
    {
        COMMON_THROWEXCEPTION( "calltree.ct not open" )
    }

    // Add a cost line for the current region 

    uint64_t vals[2] = { 10, 10 };
    int scl = 0;

    CT.add( region_id, -1, scl, vals );

    theCallStack.pop_back();
}

void CallTree::finish()
{
    if ( outfile != NULL )
    {
        fprintf( outfile, "# closed by finish\n" );
        fclose( outfile );
    }
 
    outfile = NULL;
}

static Guard myGuard;

}
