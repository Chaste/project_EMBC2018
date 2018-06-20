/*

Copyright (c) 2005-2014, University of Oxford.
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

#include "AppliedForceModifier.hpp"
#include "MeshBasedCellPopulation.hpp"

#include "Debug.hpp"
#include <math.h>


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::AppliedForceModifier()
    : AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>(),
      mResetTractionsOnCells(true),
      mTractionFile("")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::~AppliedForceModifier()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
	/*
	 * We must update CellData in SetupSolve(), otherwise it will not have been
	 * fully initialised by the time we enter the main time loop.
	 */
	if(mResetTractionsOnCells && mTractionFile !="")
	{
		LoadTractionFromFile();
		UpdateCellData(rCellPopulation);
	}
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
	assert(SPACE_DIM==3);

	// See if any edges are too long and if so divide them
	double num_cells = rCellPopulation.GetNumRealCells();

	MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_cell_population = dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(rCellPopulation));
	assert(p_cell_population != NULL);

    //If there has been any remeshing recalculate the applied force from the traction file.
    if (num_cells != rCellPopulation.GetNumRealCells())
    {
        PRINT_2_VARIABLES(num_cells,rCellPopulation.GetNumRealCells());
        UpdateCellData(rCellPopulation);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
	assert(SPACE_DIM==3); // Currently assumes that SPACE_DIM = 3

	for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
		 cell_iter != rCellPopulation.End();
		 ++cell_iter)
	{
		c_vector<double, SPACE_DIM> location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

		unsigned nearest_fluid_site = UNSIGNED_UNSET;
		double distance_to_fluid_site = DBL_MAX;
		for (unsigned fluid_site_index = 0; fluid_site_index <  mAppliedPosition.size(); fluid_site_index++)
		{
			double distance = norm_2(location - mAppliedPosition[fluid_site_index]);
			if (distance < distance_to_fluid_site)
			{
				distance_to_fluid_site = distance;
				nearest_fluid_site = fluid_site_index;
			}

		}
		assert(nearest_fluid_site != UNSIGNED_UNSET);


		// Calculate the approximate area of the voronoi region around the cell by including a third of the area
		// of all surrounding triangles. Useful for turning stresses into forces.
		double voronoi_cell_area = rCellPopulation.GetVolumeOfCell(*cell_iter);

		c_vector<double,3> force = mAppliedTractions[nearest_fluid_site]*voronoi_cell_area;
		c_vector<double,3> shear_stress = mAppliedTangentTractions[nearest_fluid_site];

		assert(fabs(force[0])<1e10);
		assert(fabs(force[1])<1e10);
		assert(fabs(force[2])<1e10);
		assert(fabs(shear_stress[0])<1e10);
		assert(fabs(shear_stress[1])<1e10);
		assert(fabs(shear_stress[2])<1e10);

		// Store the force in CellData
		cell_iter->GetCellData()->SetItem("applied_force_x", force[0]);
		cell_iter->GetCellData()->SetItem("applied_force_y", force[1]);
		cell_iter->GetCellData()->SetItem("applied_force_z", force[2]);
		cell_iter->GetCellData()->SetItem("applied_force_mag", norm_2(force));
		cell_iter->GetCellData()->SetItem("voronoi_cell_area", voronoi_cell_area);
		cell_iter->GetCellData()->SetItem("applied_shear_stress_x", shear_stress[0]);
		cell_iter->GetCellData()->SetItem("applied_shear_stress_y", shear_stress[1]);
		cell_iter->GetCellData()->SetItem("applied_shear_stress_z", shear_stress[2]);
		cell_iter->GetCellData()->SetItem("applied_shear_stress_mag", norm_2(shear_stress));
	}
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::GetResetTractionsOnCells()
{
    return mResetTractionsOnCells;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::SetResetTractionsOnCells(bool resetTractionsOnCells, std::string tractionFile)
{
	assert(!resetTractionsOnCells || tractionFile!="");
	mResetTractionsOnCells = resetTractionsOnCells;
	mTractionFile = tractionFile;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::LoadTractionFromFile()
{
	PRINT_VARIABLE(mTractionFile);
	FILE* traction_file = fopen((char*)mTractionFile.c_str(), "r");
	assert(traction_file != NULL);
	hemelb::io::writers::xdr::XdrFileReader reader(traction_file);

	// File format described in http://pauli.chem.ucl.ac.uk/trac/hemelb/wiki/ExtractionFiles

    unsigned hemelb_magic_number, extraction_magic_number, extraction_version_number;
    reader.readUnsignedInt(hemelb_magic_number);
    reader.readUnsignedInt(extraction_magic_number);
    reader.readUnsignedInt(extraction_version_number);

    assert(hemelb_magic_number == 0x686c6221);
    assert(extraction_magic_number == 0x78747204);
    assert(extraction_version_number == 4);

    double voxel_size;
    reader.readDouble(voxel_size);

    c_vector<double,3> origin;
    reader.readDouble(origin[0]);
    reader.readDouble(origin[1]);
    reader.readDouble(origin[2]);

    unsigned long long number_fluid_sites;
    reader.readUnsignedLong(number_fluid_sites);

    unsigned field_count;
    reader.readUnsignedInt(field_count);
    assert(field_count == 2); // Traction and tangetial component of traction

    unsigned header_length;
    reader.readUnsignedInt(header_length);

    // Traction field header
    std::string field_name;
    unsigned number_floats;
    double traction_offset;

    reader.readString(field_name, header_length);
    reader.readUnsignedInt(number_floats);
    assert(number_floats == 3);
    reader.readDouble(traction_offset);

    // Tangential traction field header
    double tanget_traction_offset;

    reader.readString(field_name, header_length);
    reader.readUnsignedInt(number_floats);
    assert(number_floats == 3);
    reader.readDouble(tanget_traction_offset);

    // Data section (we are reading a single timestep)
    unsigned long long timestep_num;
    reader.readUnsignedLong(timestep_num);

    mAppliedPosition.clear();
    mAppliedTractions.clear();
    mAppliedTangentTractions.clear();

    for (unsigned fluid_site_index = 0; fluid_site_index <  number_fluid_sites; fluid_site_index++)
    {
    	{
    		c_vector<unsigned,3> coords;
    		reader.readUnsignedInt(coords[0]);
    		reader.readUnsignedInt(coords[1]);
    		reader.readUnsignedInt(coords[2]);

    		mAppliedPosition.push_back(origin+voxel_size*coords);
    	}

    	{
			c_vector<float,3> traction;
			reader.readFloat(traction[0]);
			traction[0] += traction_offset;
			reader.readFloat(traction[1]);
			traction[1] += traction_offset;
			reader.readFloat(traction[2]);
			traction[2] += traction_offset;

			assert(fabs(traction[0])<1e10);
			assert(fabs(traction[1])<1e10);
			assert(fabs(traction[2])<1e10);

			mAppliedTractions.push_back(traction);
    	}

    	{
			c_vector<float,3> tangent_traction;
			reader.readFloat(tangent_traction[0]);
			tangent_traction[0] += tanget_traction_offset;
			reader.readFloat(tangent_traction[1]);
			tangent_traction[1] += tanget_traction_offset;
			reader.readFloat(tangent_traction[2]);
			tangent_traction[2] += tanget_traction_offset;

			assert(fabs(tangent_traction[0])<1e10);
			assert(fabs(tangent_traction[1])<1e10);
			assert(fabs(tangent_traction[2])<1e10);

			mAppliedTangentTractions.push_back(tangent_traction);
    	}
    }

    assert(mAppliedPosition.size() == number_fluid_sites);
    assert(mAppliedTractions.size() == number_fluid_sites);
    assert(mAppliedTangentTractions.size() == number_fluid_sites);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForceModifier<ELEMENT_DIM,SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}


// Explicit instantiation
template class AppliedForceModifier<1,1>;
template class AppliedForceModifier<1,2>;
template class AppliedForceModifier<2,2>;
template class AppliedForceModifier<1,3>;
template class AppliedForceModifier<2,3>;
template class AppliedForceModifier<3,3>;
// Serialization for Boost >= 1.36

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AppliedForceModifier)

