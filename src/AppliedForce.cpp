#include "AppliedForce.hpp"
#include "MeshBasedCellPopulation.hpp"
#include <cmath>


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AppliedForce<ELEMENT_DIM, SPACE_DIM>::AppliedForce()
    : AbstractForce<ELEMENT_DIM, SPACE_DIM>()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AppliedForce<ELEMENT_DIM, SPACE_DIM>::~AppliedForce()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForce<ELEMENT_DIM, SPACE_DIM>::AddForceContribution(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
	assert(SPACE_DIM==3); // Currently assumes that SPACE_DIM = 3

    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned node_global_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

        c_vector<double, SPACE_DIM> applied_force;

        // Retrieve the force from CellData and add to rForces
        applied_force[0]=cell_iter->GetCellData()->GetItem("applied_force_x");
        applied_force[1]=cell_iter->GetCellData()->GetItem("applied_force_y");
        applied_force[2]=cell_iter->GetCellData()->GetItem("applied_force_z");
        
        assert(fabs(applied_force[0])<1e10);
        assert(fabs(applied_force[1])<1e10);
        assert(fabs(applied_force[2])<1e10);

        rCellPopulation.GetNode(node_global_index)->AddAppliedForceContribution(applied_force);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AppliedForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // No parameters to include

    // Call method on direct parent class
    AbstractForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////
template class AppliedForce<2,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AppliedForce)
