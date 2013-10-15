#include <unistd.h>
#include <cstdlib>
#include <cstring>
#include <stdarg.h>
#include <iostream>
#include <vector>

#define PNG_DEBUG 3

#include <lama.hpp>
#include <png.h>

#include <lama/storage/CSRStorage.hpp>
#include <lama/HostReadAccess.hpp>

class Bitmap
{
public:

    /** Generate an empty bitmap of a certain size.
     *
     *  @param[in] w width of the bitmap
     *  @param[in] h height of the bitmap
     *  @param[in] s scaling, e.g, 1, 2, .. 10
     *
     *  One pixel corresponds will be presented by s x s pixels in the final written bitmap.
     */

    Bitmap( const int w, const int h, const int s )
    {
        LAMA_ASSERT( s > 0, "scale factor s = " << s << " must be positive" )

        scale  = s;
        width  = w;
        height = h;

        bit_depth  = 8;
        color_type = PNG_COLOR_TYPE_RGBA;

        row_pointers = new png_bytep[ height * scale ];

        Color background( 220, 220, 220 );

        for ( int y = 0; y < height * scale; ++y )
        {
            row_pointers[y] = new png_byte[ 4 * width * scale ];

            for ( int x = 0; x < width * scale; ++x )
            {
               row_pointers[y][4 * x ] = background.r;
               row_pointers[y][4 * x + 1 ] = background.g;
               row_pointers[y][4 * x + 2 ] = background.b;
               row_pointers[y][4 * x + 3 ] = 255;
            }
        }

        setMinMax( 0, 1 );
    }

    ~Bitmap()
    {
        for ( int y = 0; y < height * scale; y++ )
        {
            delete [] row_pointers[y];
        }

        delete [] row_pointers;
    }

    void setPixel( int x, int y, int r, int g, int b )
    {
        row_pointers[x][4 * y ] = r;
        row_pointers[x][4 * y + 1 ] = g;
        row_pointers[x][4 * y + 2 ] = b;
    }

    void set( int x, int y, double val )
    {
        if ( x < 0 || x >= height )
        {
            std::cerr << "x = " << x << " out of range: 0 .. " << height << std::endl;
            exit( 1 );
        }

        if ( y < 0 || y >= width )
        {
            std::cerr << "y = " << y << " out of range: 0 .. " << width << std::endl;
            exit( 1 );
        }

        int r, g, b;

        getColor( r, g, b, val );

        for ( int i = 0; i < scale; ++i )
        {
            for ( int j = 0; j < scale; ++j )
            { 
                setPixel( scale * x + i, scale * y + j, r, g, b );
            }
        }
    }

    template<typename ValueType>
    void drawCSR( const int nRows, const int nCols, const int ia[], const int ja[], const ValueType values[] )
    {
        ValueType minval =  10000.0;
        ValueType maxval = -100000.0;

        for ( int k = 0; k < ia[nRows]; ++k )
        {
            ValueType v = values[k];
            if ( v < minval ) minval = v;
            if ( v > maxval ) maxval = v;
        }

        double multRow = double( height ) / double( nRows );
        double multCol = double( width ) / double( nCols );

        setMinMax( minval, maxval);
 
        for ( int i = 0; i < nRows; ++i )
        {
            for ( int j = ia[i]; j < ia[i+1]; ++j )
            {
                set( static_cast<int>( i * multRow ), static_cast<int>( ja[j] * multCol ), values[j] );
            }
        }

        // diagonals drawn at the end

        for ( int i = 0; i < nRows; ++i )
        {
            for ( int j = ia[i]; j < ia[i+1]; ++j )
            {
                if ( ja[j] == i )
                {
                    set( static_cast<int>( i * multRow ), static_cast<int>( ja[j] * multCol ), values[j] );
                }
            }
        }
    }

    void getColor( int& r, int& g, int& b, double val )
    {
        if ( palette.size() == 0 )
        {
            std::cerr << "No colors defined" << std::endl;
            exit(1);
        }
        else if ( palette.size() == 1 )
        {
            r = palette[0].r;
            g = palette[0].g;
            b = palette[0].b;
        }
        else 
        {
            double sval = val;
            if ( sval < minval) sval = minval;
            if ( sval > maxval) sval = maxval;

            sval = ( sval - minval ) / ( maxval - minval ); // normed between 0.0 and 1.0

            r = static_cast<int>( palette[0].r + sval * ( palette[1].r - palette[0].r ) );
            g = static_cast<int>( palette[0].g + sval * ( palette[1].g - palette[0].g ) );
            b = static_cast<int>( palette[0].b + sval * ( palette[1].b - palette[0].b ) );
        }
    }

    /** Write the picture data as a png file. */

    void write_png_file( const char* file_name )
    {
        /* create file */
        FILE* fp = fopen( file_name, "wb" );

        if ( !fp )
        {
            std::cerr << "[write_png_file] File " << file_name << " could not be opened for writing" << std::endl;
            exit( 1 );
        }

        /* initialize stuff */
        png_ptr = png_create_write_struct( PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );

        if ( !png_ptr )
        {
            std::cerr << "[write_png_file] png_create_write_struct failed" << std::endl;
            exit( 1 );
        }

        info_ptr = png_create_info_struct( png_ptr );

        if ( !info_ptr )
        {
            std::cerr << "[write_png_file] png_create_info_struct failed" << std::endl;
            exit( 1 );
        }

        if ( setjmp( png_jmpbuf( png_ptr ) ) )
        {
            std::cerr << "[write_png_file] Error during init_io" << std::endl;
        }

        png_init_io( png_ptr, fp );

        /* write header */

        if ( setjmp( png_jmpbuf( png_ptr ) ) )
        {
            std::cerr << "[write_png_file] Error during writing header" << std::endl;
        }

        png_set_IHDR( png_ptr, info_ptr, width * scale, height * scale,
                      bit_depth, color_type, PNG_INTERLACE_NONE,
                      PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE );

        png_write_info( png_ptr, info_ptr );

        /* write bytes */
        if ( setjmp( png_jmpbuf( png_ptr ) ) )
        {
            std::cerr << "[write_png_file] Error during writing bytes" << std::endl;
            exit( 1 );
        }

        png_write_image( png_ptr, row_pointers );

        /* end write */

        if ( setjmp( png_jmpbuf( png_ptr ) ) )
        {
            std::cerr << "[write_png_file] Error during end of write" << std::endl;
            exit( 1 );
        }

        png_write_end( png_ptr, NULL );
        fclose( fp );
    }

    /** Set color for drawing matrix entries. */

    void setColor( const int r, const int g, const int b )
    {
        Color color( r, g, b );
        palette.push_back( color );
    }

private:

    void setMinMax( double low, double high )
    {
        minval = low;
        maxval = high;
    }

    int width;
    int height;

    int scale;

    png_byte color_type;
    png_byte bit_depth;

    png_structp png_ptr;
    png_infop info_ptr;
    png_bytep* row_pointers;

    double minval;
    double maxval;

    struct Color
    {
        int r, g, b;
 
        Color( const int red, const int green, const int blue )
        {
            r = red;
            g = green;
            b = blue;
        }
    };
   
    std::vector<Color> palette;
};

using namespace lama;

int main( int argc, char** argv )
{
    CSRStorage<double> matrix;

    if ( argc < 2 )
    {
        std::cerr << "Missing filename for input matrix" << std::endl;
        std::cerr << "spy matrix_filename [ width [ height [ scale ] ] ]" << std::endl;
        exit(1);
    }

    const char* filename = argv[1];

    int nRows = 800;

    if ( argc > 2 )
    {
        sscanf( argv[2], "%d",  &nRows );
    }

    int nColumns = nRows;

    if ( argc > 3 )
    {
        sscanf( argv[3], "%d",  &nColumns );
    }

    int nZoom = 1;

    if ( argc > 4 )
    {
        sscanf( argv[4], "%d",  &nZoom );
    }

    matrix.readFromFile( filename );

    const LAMAArray<IndexType>& ia = matrix.getIA();
    const LAMAArray<IndexType>& ja = matrix.getJA();
    const LAMAArray<double>& values = matrix.getValues();

    HostReadAccess<IndexType> csrIA( ia );
    HostReadAccess<IndexType> csrJA( ja );
    HostReadAccess<double> csrValues( values );

    std::cout << "Write png of size " << nRows << " x " << nColumns << ", zoom = " << nZoom << std::endl;

    Bitmap pic( nRows, nColumns, nZoom );

    pic.setColor( 240, 120, 0 );  // color for smallest value
    // pic.setColor( 0, 0, 255 );    // color for largetst value

    pic.drawCSR( matrix.getNumRows(), matrix.getNumColumns(), csrIA.get(), csrJA.get(), csrValues.get() );

    const std::string out_filename = "lama.png";

    pic.write_png_file( out_filename.c_str() );

    std::cout << "png files has been written as " << out_filename << std::endl;
}
