#ifndef APPLIEDFORCE_HPP_
#define APPLIEDFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"

// Needed here to avoid serialization errors (on Boost<1.37)
#include "CellLabel.hpp"

/**
 * A force class to be used with an AppliedForceOffLatticeSimulation to impose an externally defined force.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AppliedForce  : public AbstractForce<ELEMENT_DIM, SPACE_DIM>
{
friend class TestForces;

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /**
     * Constructor.
     */
    AppliedForce();

    /**
     * Destructor.
     */
    ~AppliedForce();

    /**
     * Overridden AddForceContribution() method. Which uses the applied force stored in CellData.
     *
     * @param rCellPopulation reference to the cell population
     *
     */
    void AddForceContribution(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AppliedForce)

#endif /*APPLIEDFORCE_HPP_*/
