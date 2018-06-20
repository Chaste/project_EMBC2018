#ifndef FIXEDREGIONBOUNDARYCONDITION_HPP_
#define FIXEDREGIONBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A plane cell population boundary condition class, which stops nodes moving through
 * a specified plane in the domain. Although the name of this class suggests it is
 * specific to 3D, it is actually also implemented for 1D and 2D, for which it is
 * really a 'point' and 'line' boundary condition respectively.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class FixedRegionBoundaryCondition : public AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM>
{
private:

    /**
     * A point on the boundary plane.
     */
    c_vector<double, SPACE_DIM> mPointOnPlane;

    /**
     * The outward-facing unit normal vector to the boundary plane.
     */
    c_vector<double, SPACE_DIM> mNormalToPlane;
    
    /**
     * The maximal radius of influence of the boundary condition
     */
    double mRadius;
    
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> >(*this);
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param point a point on the boundary plane
     * @param normal the outward-facing unit normal vector to the boundary plane
     * @param radius the maximal influence of the boundary condition from the point
     */
    FixedRegionBoundaryCondition(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>* pCellPopulation,
                           c_vector<double, SPACE_DIM> point,
                           c_vector<double, SPACE_DIM> normal,
                           double radius);

    /**
     * @return mPointOnPlane.
     */
    const c_vector<double, SPACE_DIM>& rGetPointOnPlane() const;

    /**
     * @return mNormalToPlane.
     */
    const c_vector<double, SPACE_DIM>& rGetNormalToPlane() const;
    
    /**
     * @return mRadius.
     */
    const double& rGetRadius() const;
    
    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> >& rOldLocations);

    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition();

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(FixedRegionBoundaryCondition)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a FixedRegionBoundaryCondition.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const FixedRegionBoundaryCondition<ELEMENT_DIM,SPACE_DIM>* t, unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
    // Archive c_vectors one component at a time
    c_vector<double, SPACE_DIM> point = t->rGetPointOnPlane();
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        ar << point[i];
    }
    c_vector<double, SPACE_DIM> normal = t->rGetNormalToPlane();
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        ar << normal[i];
    }
    
    ar << t->rGetRadius();
}

/**
 * De-serialize constructor parameters and initialize a FixedRegionBoundaryCondition.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, FixedRegionBoundaryCondition<ELEMENT_DIM,SPACE_DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_cell_population;
    ar >> p_cell_population;

    // Archive c_vectors one component at a time
    c_vector<double, SPACE_DIM> point;
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        ar >> point[i];
    }
    c_vector<double, SPACE_DIM> normal;
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        ar >> normal[i];
    }
    double radius;
    ar >> radius;

    // Invoke inplace constructor to initialise instance
    ::new(t)FixedRegionBoundaryCondition<ELEMENT_DIM,SPACE_DIM>(p_cell_population, point, normal, radius);
}
}
} // namespace ...

#endif /*FIXEDREGIONBOUNDARYCONDITION_HPP_*/
