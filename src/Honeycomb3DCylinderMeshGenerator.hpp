/*

Copyright (c) 2005-2013, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef HONEYCOMB3DCYLINDRICALMESHGENERATOR_HPP_
#define HONEYCOMB3DCYLINDRICALMESHGENERATOR_HPP_

#include <cmath>
#include <vector>

#include "MutableMesh.hpp"

/**
 * Honeycomb mesh generator that creates a 2D honeycomb surface mesh (with equal distance
 * between nodes) for use in Vascular remodelling simulations.
 *
 * NOTE: the user should delete the mesh after use to manage memory.
 */
class Honeycomb3DCylinderMeshGenerator
{
protected:

    /** A pointer to the mesh this class creates */
    MutableMesh<2,3>* mpMesh;

    /** The mesh is generated by writing out a series of nodes and reading them in from this file*/
    std::string mMeshFilename;

public:

    /**
     * Default constructor.
     *
     * @param numNodesAround  The number of cells you want around the cylinder
     * @param numNodesAlongLength  The number of cells you want sides of the domain
     * @param radius radius of cylinder
     * @param scaleFactor length of cylinder
     */
    Honeycomb3DCylinderMeshGenerator(unsigned numNodesAround, unsigned numNodesAlongLength, double radius, double length);

    /**
     * Null constructor for derived classes to call.
     */
    Honeycomb3DCylinderMeshGenerator()
    {
    }

    /**
     * Destructor - deletes the mesh object and pointer
     */
    virtual ~Honeycomb3DCylinderMeshGenerator();

    /**
     * @return a the cylinder mesh
     */
    virtual MutableMesh<2,3>* GetMesh();
};

#endif /*HONEYCOMB3DCYLINDRICALMESHGENERATOR_HPP_*/
