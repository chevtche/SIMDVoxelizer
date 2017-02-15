/* Copyright (c) 2015-2016, EPFL/Blue Brain Project
 * All rights reserved. Do not distribute without permission.
 * Responsible Author: Cyrille Favreau <cyrille.favreau@epfl.ch>
 *
 * This file is part of SIMDVoxelizer <https://github.com/favreau/SIMDVoxelizer>
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License version 3.0 as published
 * by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include "SIMDSparseVoxelizer_ispc.h"
#include "Octree.h"
#include <string.h>
#include <cmath>
#include <fstream>
#include <vector>
#include <map>
#include <limits>
#include <iostream>
#include <sstream>
#include <ctime>
#include <cstdlib>

using namespace ispc;

int span = 32;
float voxelSize = 2.f;
std::string inputFile;
std::string outputFile;

int main( int argc, char* argv[] )
{
    if( argc != 5 )
    {
        std::cerr << "usage: SIMDVoxelizer <voxel_size> <cutoff_distance> "
                  << "<input_file> <output_file>" << std::endl;
        exit(1);
    }

    voxelSize = atof(argv[1]);
    span = atoi(argv[2]);
    inputFile = argv[3];
    outputFile = argv[4];
    std::vector<float> events;
    std::ifstream file(inputFile.c_str(), std::ios::in | std::ios::binary);
    if( file.is_open( ))
    {
        while( !file.eof( ))
        {
            float v;
            file.read( (char *)&v, sizeof(float));
            events.push_back(v);
        }
        file.close();
    }

    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "SIMDVoxelizer" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "Voxel size        : " << voxelSize << std::endl;
    std::cout << "Span              : " << span << std::endl;
    std::cout << "Input file        : " << inputFile << std::endl;
    std::cout << "Output file       : " << outputFile << std::endl;
    std::cout << "--------------------------------------------" << std::endl;

    Octree morphoOctree( events, voxelSize );

    uint64_t volumeSize = morphoOctree.getVolumeSize();
    glm::uvec3 volumeDim = morphoOctree.getVolumeDim();

    float* volume = new float[ volumeSize ];

    std::cout<< "Volume dims: " << volumeDim.x << " " << volumeDim.y << " " << volumeDim.z << " " << volumeSize << std::endl;

    uint32_t zLenght = 32;
    uint32_t nPasses = std::ceil( volumeDim.z / zLenght );

    for( uint32_t zOffset = 0; zOffset < volumeDim.z; zOffset += zLenght )
    {
        std::cout << "z: " << zOffset << std::endl;
        SIMDSparseVoxelizer_ispc( zOffset, zOffset + zLenght, span, voxelSize,
                                  morphoOctree.getOctreeSize(), volumeDim.x,
                                  volumeDim.y, volumeDim.z,
                                  morphoOctree.getFlatIndexes(),
                                  morphoOctree.getFlatData(), volume );
    }

    float minValue = std::numeric_limits<float>::max();
    float maxValue = -std::numeric_limits<float>::max();
    for( int i = 0; i < volumeSize; ++i )
    {
        if( volume[i] != 0.f )
        {
            minValue = std::min( minValue, volume[i] );
            maxValue = std::max( maxValue, volume[i] );
        }
    }

    //float a = 255.f / ( maxValue - minValue );
    float a = 9.f / ( maxValue - minValue );
    std::cout << "Normalization [" << minValue << " - " << maxValue << "] " << a << std::endl;

    char* volumeAsChar = new char[ volumeSize ];
    for( int i = 0; i < volumeSize; ++i )
    {
        //float normalizedValue = (volume[i] - minValue) * a;
        float normalizedValue = 255.0f * std::log((volume[i] - minValue) * a + 1.0);
        volumeAsChar[i] = (uint8_t)normalizedValue;
    }
    std::cout << std::endl;
    std::cout << "--------------------------------------------" << std::endl;

    std::ofstream volumeFile(
        outputFile.c_str(), std::ios::out | std::ios::binary);
    volumeFile.write( (char*)&volumeAsChar[0], sizeof(char) * volumeSize );
    volumeFile.close();
    delete [] volume;
    delete [] volumeAsChar;

    return 0;
}
