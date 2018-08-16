/**
 * @file solver/examples/FDSimulation/Configuration.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief ToDo: Missing description in ./solver/examples/FDSimulation/Configuration.hpp
 * @author Lauretta Schubert
 * @date 27.09.2016
 */
#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>

using namespace scai;

template<typename ValueType>
class Configuration
{
public:
    Configuration( std::string filename )
    {
        // read all lines in file

        std::string line;
        std::map<std::string, std::string> map;
        std::ifstream input( filename.c_str() );

        while ( std::getline( input, line ) )
        {
            size_t lineEnd = line.size();
            std::string::size_type commentPos1 = line.find_first_of( "/", 0 );

            if ( std::string::npos != commentPos1 )
            {
                std::string::size_type commentPos2 = line.find_first_of( "/", commentPos1 );

                if ( std::string::npos != commentPos2 )
                {
                    if ( commentPos1 == 0 )
                    {
                        continue;
                    }

                    lineEnd = commentPos1;
                }
            }

            std::string::size_type equalPos = line.find_first_of( "=", 0 );

            if ( std::string::npos != equalPos )
            {
                // tokenize it  name = val
                std::string name = line.substr( 0, equalPos );
                size_t len = lineEnd - ( equalPos + 1 );
                std::string val  = line.substr( equalPos + 1, len );
                map.insert( std::pair<std::string, std::string>( name, val ) );
            }
        }

        input.close();

        // check map and assign all values with right "cast" to members
        size_t nArgs = 16;

        if ( map.size() != nArgs )
        {
            std::cout << filename << " does not include a valid configutation with " << nArgs << " arguments." << std::endl;
        }

        std::istringstream( map[ "NZ" ] ) >> NZ; // IndexType
        std::istringstream( map[ "NX" ] ) >> NX; // IndexType
        std::istringstream( map[ "NY" ] ) >> NY; // IndexType

        std::istringstream( map[ "DH" ] ) >> DH; // IndexType

        std::istringstream( map[ "DT" ] ) >> DT; // ValueType
        std::istringstream( map[ "T" ] ) >> T;  // ValueType

        std::istringstream( map[ "velocity" ] ) >> velocity; // ValueType
        std::istringstream( map[ "rho" ] ) >> rho; // ValueType

        std::istringstream( map[ "fc" ] ) >> fc; // IndexType
        std::istringstream( map[ "amp" ] ) >> amp; // IndexType

        std::istringstream( map[ "source_z" ] ) >> source_z; // IndexType
        std::istringstream( map[ "source_x" ] ) >> source_x; // IndexType
        std::istringstream( map[ "source_y" ] ) >> source_y; // IndexType

        std::istringstream( map[ "seismogram_z" ] ) >> seismogram_z; // IndexType
        std::istringstream( map[ "seismogram_x" ] ) >> seismogram_x; // IndexType
        std::istringstream( map[ "seismogram_y" ] ) >> seismogram_y; // IndexType

        // calculate other parameters

        N = NZ * NX * NY;

        M = velocity * velocity * rho; // P-wave modulus

        NT = static_cast<IndexType>( ( T / DT ) + 0.5 ); // MATLAB round(T/DT)

        v_factor = DT / DH / rho;
        p_factor = DT * M;

        source_index     = index( source_x,     source_y,     source_z,     NX, NY, NZ );
        seismogram_index = index( seismogram_x, seismogram_y, seismogram_z, NX, NY, NZ );

    }

    ~Configuration() {}

    /*
     *  routine for printing out the configuration
     */
    void print()
    {
        IndexType velocity_max = static_cast<IndexType>( velocity ); // TODO: is velocity a vector?
        double courant = velocity_max * DT / DH;
        IndexType wavelength = velocity_max / fc;

        std::cout << "Configuration:" << std::endl;
        std::cout << "Criteriums:" << std::endl;
        std::cout << "    Wavelength: " << wavelength << " m" << std::endl;
        std::cout << "    Gridpoints per Wavelength: " << wavelength / DH << std::endl;
        std::cout << "    Courant-number: " << courant << std::endl;

        if ( courant >= 0.8 )
        {
            std::cout << "Simulation will be UNSTABLE" << std::endl;
            std::cout << "Choose smaller DT, eg.: " << DH * 0.3 / velocity_max << std::endl;
            exit( 0 );
        }

        std::cout << "Modelling-domain:" << std::endl;
        std::cout << "    Z: " << DH* NZ << " m (Depth)" << std::endl;
        std::cout << "    X: " << DH* NX << " m (Horizontal)" << std::endl;
        std::cout << "    Y: " << DH* NY << " m (Horizontal)" << std::endl;
        std::cout << "Material:" << std::endl;
        std::cout << "    Velocity:" << velocity << " m/s" << std::endl;
        std::cout << std::endl;
    }

    // for debug
    void printAllRaw()
    {
        std::cout << "NZ=" << getNZ() << std::endl;
        std::cout << "NX=" << getNX() << std::endl;
        std::cout << "NY=" << getNY() << std::endl;
        std::cout << "DH=" << getDH() << std::endl;
        std::cout << "DT=" << getDT() << std::endl;
        std::cout << "T=" << getT() << std::endl;
        std::cout << "velocity=" << getVelocity() << std::endl;
        std::cout << "rho=" << getRho() << std::endl;
        std::cout << "fc=" << getFC() << std::endl;
        std::cout << "amp=" << getAMP() << std::endl;
        std::cout << "source_z=" << getSourceZ() << std::endl;
        std::cout << "source_x=" << getSourceX() << std::endl;
        std::cout << "source_y=" << getSourceY() << std::endl;
        std::cout << "seismogram_z=" << getSeismogramZ() << std::endl;
        std::cout << "seismogram_x=" << getSeismogramX() << std::endl;
        std::cout << "seismogram_y=" << getSeismogramX() << std::endl;
        std::cout << "N=" << getN() << std::endl;
        std::cout << "M=" << getM() << std::endl;
        std::cout << "NT=" << getNT() << std::endl;
        std::cout << "v_factor=" << getVfactor() << std::endl;
        std::cout << "p_factor=" << getPfactor() << std::endl;
        std::cout << "source_index=" << getSourceIndex() << std::endl;
        std::cout << "seismogram_index=" << getSeismogramIndex() << std::endl;
    }

    IndexType getNZ()
    {
        return NZ;
    }
    IndexType getNX()
    {
        return NX;
    }
    IndexType getNY()
    {
        return NY;
    }

    IndexType getDH()
    {
        return DH;
    }

    ValueType getDT()
    {
        return DT;
    }
    ValueType getT()
    {
        return T;
    }

    ValueType getVelocity()
    {
        return velocity;
    }
    ValueType getRho()
    {
        return rho;
    }

    IndexType getFC()
    {
        return fc;
    }
    IndexType getAMP()
    {
        return amp;
    }

    IndexType getSourceZ()
    {
        return source_z;
    }
    IndexType getSourceX()
    {
        return source_x;
    }
    IndexType getSourceY()
    {
        return source_y;
    }

    IndexType getSeismogramZ()
    {
        return seismogram_z;
    }
    IndexType getSeismogramX()
    {
        return seismogram_x;
    }
    IndexType getSeismogramY()
    {
        return seismogram_y;
    }

    IndexType getN()
    {
        return N;
    }

    ValueType getM()
    {
        return M;
    }

    IndexType getNT()
    {
        return NT;
    }

    ValueType& getVfactor()
    {
        return v_factor;
    }
    ValueType& getPfactor()
    {
        return p_factor;
    }

    IndexType getSourceIndex()
    {
        return source_index;
    }
    IndexType getSeismogramIndex()
    {
        return seismogram_index;
    }

private:
    /*
     *  routine for calculating the 3D index position
     *  order of directions in 3d space: z, x, y
     */
    IndexType index( IndexType x, IndexType y, IndexType z, IndexType NX, IndexType NY, IndexType NZ )
    {
        SCAI_REGION( "Index_calculation" )

        if ( z > NZ || x > NX || y > NY || z < 1 || x < 1 || y < 1 )
        {
            COMMON_THROWEXCEPTION ( "Could not map from coordinate to indize!" )
            return invalidIndex;
        }
        else
        {
            return ( ( z - 1 ) + ( x - 1 ) * NZ + ( y - 1 ) * NZ * NX );
        }
    }

    /* read parameters */

    // define spatial sampling: number of grid points in direction
    IndexType NZ; // depth
    IndexType NX; // horizontal 1
    IndexType NY; // horizontal 2

    // define distance between two grid points in meter
    IndexType DH;

    // define temporal sampling
    ValueType DT;  // temporal sampling in seconds
    ValueType T;   // total simulation time

    // define material parameter
    ValueType velocity; // Density in kilo gramms per cubic meter
    ValueType rho;      // P-wave velocity in meter per seconds

    // define source wavelet
    IndexType fc; // Center frequency of ricker wavelet
    IndexType amp; // Amplitude of source signal

    // source position in grid points
    IndexType source_z;
    IndexType source_x;
    IndexType source_y;

    // seismogram position in grid points
    IndexType seismogram_z;
    IndexType seismogram_x;
    IndexType seismogram_y;

    /* calculated parameters */

    IndexType N;

    ValueType M; // P-wave modulus

    IndexType NT;

    ValueType v_factor;
    ValueType p_factor;

    IndexType source_index;
    IndexType seismogram_index;
};
