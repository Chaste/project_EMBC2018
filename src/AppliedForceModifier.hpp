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

#ifndef APPLIEDFORCEMODIFIER_HPP_
#define APPLIEDFORCEMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"
#include "PetscTools.hpp"
#include "XdrFileReader.hpp"
#include <map>

/**
 * A modifier class which at each simulation time step
 *
 * To be used in conjunction with the AppliedForce Class
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AppliedForceModifier : public AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>
{
private:

	/**
	 * 	Store the loaded traction data; this is the applied traction at a given point.
	 */
	std::vector < c_vector<double,3> > mAppliedTractions;

	/**
	 * 	Store the loaded tangential component of the traction data; this is the applied tangential traction at a given point.
	 */
	std::vector < c_vector<double,3> > mAppliedTangentTractions;

	/**
	 * 	Store the loaded traction data; this is the location that the traction is defined at.
	 */
	std::vector < c_vector<double,3> > mAppliedPosition;

	bool mResetTractionsOnCells;

	std::string mTractionFile;

	/** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM> >(*this);
        archive & mAppliedTractions;
        archive & mAppliedTangentTractions;
        archive & mAppliedPosition;
        archive & mResetTractionsOnCells;
        archive & mTractionFile;
    }

public:

    /**
	 * Default constructor.
	 */
    AppliedForceModifier();

    /**
     * Destructor.
     */
    virtual ~AppliedForceModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory);

    /**
	 * Return the force applied to a point in space.
	 *
	 * @param location the location in space where the force is to be calculated
	 */
	c_vector<double, SPACE_DIM> GetAppliedForce(c_vector<double, SPACE_DIM> location);

	/**
	 * Helper method to load the traction data from file.
	 */
	void LoadTractionFromFile();

	/**
	 * @return mResetTractionsOnCells
	 */
	bool GetResetTractionsOnCells();

	/**
	 * Set mResetTractionsOnCells.
	 *
	 * @param resetTractionsOnCells the new value of mResetTractionsOnCells
	 */
	void SetResetTractionsOnCells(bool resetTractionsOnCells, std::string tractionFile="");

	double GetEdgeDivisionThreshold();

	void SetEdgeDivisionThreshold(double edgeDivisionThreshold);

    /**
     * Helper method to store the applied tractions in CellData.
     */
    void UpdateCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AppliedForceModifier);

#endif /*APPLIEDFORCEMODIFIER_HPP_*/
