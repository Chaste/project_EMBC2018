
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

#include "Honeycomb3DCylinderMeshGenerator.hpp"

#include <boost/foreach.hpp>
#include "TrianglesMeshReader.hpp"
#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
#include "MathsCustomFunctions.hpp"
#include "ChasteSyscalls.hpp"
#include "Debug.hpp"

Honeycomb3DCylinderMeshGenerator::Honeycomb3DCylinderMeshGenerator(unsigned numNodesAround, unsigned numNodesAlongLength, double radius, double length)
  : mpMesh(NULL),
    mMeshFilename("mesh")
{
    // The code below won't work in parallel
    assert(PetscTools::IsSequential());

    // Get a unique temporary foldername
    std::stringstream pid;
    pid << getpid();
    OutputFileHandler output_file_handler("cylinder_temporary_honeycomb_mesh_" + pid.str());
    std::string output_dir = output_file_handler.GetOutputDirectoryFullPath();

    double angle_spacing = 2.0*M_PI / (double)numNodesAround;
    double length_spacing = length / ((double)numNodesAlongLength -1);

    unsigned num_nodes             = numNodesAround*numNodesAlongLength;
	unsigned num_face_around       = numNodesAround;
	unsigned num_face_along_length = numNodesAlongLength-1;
	unsigned num_face              = 2*num_face_around*num_face_along_length;
	unsigned num_edges             = 3*num_face_around*num_face_along_length + num_face_around + num_face_along_length;

	double z0 = 0;
	double theta0 = 0.0;

	// Write node file
	out_stream p_node_file = output_file_handler.OpenOutputFile(mMeshFilename+".node");
	(*p_node_file) << std::scientific;
	//(*p_node_file) << std::setprecision(20);
	(*p_node_file) << num_nodes << "\t3\t0\t1" << std::endl;

	unsigned node = 0;
	for (unsigned i=0; i<numNodesAlongLength; i++)
	{
		for (unsigned j=0; j<numNodesAround; j++)
		{

			unsigned boundary = 0;
			if ((i==0) || (i==numNodesAlongLength-1))
			{
				boundary = 1;
			}

			double theta = theta0 + angle_spacing*((double)j + 0.25*(1.0+ SmallPow(-1.0,i+1)));
			double x = radius * cos (theta);
			double y = radius * sin (theta);
			double z = z0 + length_spacing* ((double)i);

			(*p_node_file) << node++ << "\t" << x << "\t" << y << "\t" << z << "\t" << boundary << std::endl;
		}
	}
	p_node_file->close();

	// Write faceent file and edge file
	out_stream p_face_file = output_file_handler.OpenOutputFile(mMeshFilename+".face");
	(*p_face_file) << std::scientific;

	out_stream p_edge_file = output_file_handler.OpenOutputFile(mMeshFilename+".edge");
	(*p_node_file) << std::scientific;

	(*p_face_file) << num_face << "\t0" << std::endl;
	(*p_edge_file) << num_edges << "\t1" << std::endl;

	unsigned face = 0;
	unsigned edge = 0;
	for (unsigned i=0; i<num_face_along_length; i++)
	{
		for (unsigned j=0; j < num_face_around; j++)
		{
			unsigned node0 =     i*numNodesAround + j;
			unsigned node1 =     i*numNodesAround + j+1;
			unsigned node2 = (i+1)*numNodesAround + j;

			if (j == num_face_around -1)
			{
				node1 = node1 - numNodesAround;

			}
			if (i%2 != 0)
			{
				node2 = node2 + 1;
				if (j == num_face_around -1)
				{
					node2 = node2 - numNodesAround;

				}
			}

			(*p_face_file) << face++ << "\t" << node0 << "\t" << node1 << "\t" << node2 << std::endl;

			unsigned horizontal_edge_is_boundary_edge = 0;
			unsigned vertical_edge_is_boundary_edge = 0;
			if (i==0)
			{
				horizontal_edge_is_boundary_edge = 1;
			}

			(*p_edge_file) << edge++ << "\t" << node0 << "\t" << node1 << "\t" << horizontal_edge_is_boundary_edge << std::endl;
			(*p_edge_file) << edge++ << "\t" << node1 << "\t" << node2 << "\t" << 0 << std::endl;
			(*p_edge_file) << edge++ << "\t" << node2 << "\t" << node0 << "\t" << vertical_edge_is_boundary_edge << std::endl;

			node0 = i*numNodesAround + j + 1;
			node1 = (i+1)*numNodesAround + j+1;
			node2 = (i+1)*numNodesAround + j;

			if (i%2 != 0)
			{
				node0 = node0 - 1;
				if (j == num_face_around -1)
				{
					node1 = node1 - numNodesAround;
				}
			}
			else
			{
				if (j == num_face_around -1)
				{
					node0 = node0 - numNodesAround;
					node1 = node1 - numNodesAround;
				}
			}

			(*p_face_file) << face++ << "\t" << node0 << "\t" << node1 << "\t" << node2 << std::endl;
		}
	}

	for (unsigned i=0; i<num_face_along_length; i++)
	{
		unsigned node0, node1;

		if (i%2==0)
		{
			 node0 = (i+1)*numNodesAround - 1;
			 node1 = (i+2)*numNodesAround - 1;
		}
		else
		{
			node0 = (i+1)*numNodesAround;
			node1 = (i)*numNodesAround;
		}
		(*p_edge_file) << edge++ << "\t" << node0 << "\t" << node1 << "\t" << 1 << std::endl;
	}

	for (unsigned j=0; j<num_face_around; j++)
	{
		unsigned node0 = numNodesAround*(numNodesAlongLength-1) + j;
		unsigned node1 = numNodesAround*(numNodesAlongLength-1) + j+1;

		if( j== num_face_around -1 )
		{
			node1 = node1 -numNodesAround;
		}
		(*p_edge_file) << edge++ << "\t" << node1 << "\t" << node0 << "\t" << 1 << std::endl;
	}

	p_face_file->close();
	p_edge_file->close();

	// Having written the mesh to file, now construct it using TrianglesMeshReader.
	// Nested scope so the reader closes files before we delete them below.
	{
		TrianglesMeshReader<2,3> mesh_reader(output_file_handler.GetOutputDirectoryFullPath() + mMeshFilename);
		mpMesh = new MutableMesh<2,3>();
		mpMesh->ConstructFromMeshReader(mesh_reader);
	}

//	// Make the mesh cylindrical (we use Triangle library mode inside this ReMesh call)
//	mpMesh->ReMesh();
//
//	// Delete the temporary files
//	output_file_handler.FindFile("").Remove();

	// The original files have been deleted, it is better if the mesh object forgets about them
	mpMesh->SetMeshHasChangedSinceLoading();
}

Honeycomb3DCylinderMeshGenerator::~Honeycomb3DCylinderMeshGenerator()
{
    delete mpMesh;
}

MutableMesh<2,3>* Honeycomb3DCylinderMeshGenerator::GetMesh()
{
    return mpMesh;
}
