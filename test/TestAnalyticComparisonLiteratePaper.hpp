#ifndef TESTCYLINDERVALIDATION_HPP_
#define TESTCYLINDERVALIDATION_HPP_

#include <cxxtest/TestSuite.h>

#include <cstdio>
#include <ctime>
#include <cmath>

#include "OffLatticeSimulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "SmartPointers.hpp"
#include "CellLabel.hpp"
#include "VtkMeshWriter.hpp"
#include "Debug.hpp"
#include "CommandLineArguments.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "CellIdWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "PetscSetupAndFinalize.hpp"

static const double M_TIME_FOR_SIMULATION = 50.0; //50

class RadialForce : public AbstractForce<2,3>
{
private:

    double mStrength;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<2,3> >(*this);
        archive & mStrength;
    }

public:
    RadialForce(double strength=1.0)
        : AbstractForce<2,3>(),
          mStrength(strength)
    {
        assert(mStrength > 0.0);
    }

    void AddForceContribution(AbstractCellPopulation<2,3>& rCellPopulation)
    {
        // Helper variables
        MeshBasedCellPopulation<2,3>* p_cell_population = static_cast<MeshBasedCellPopulation<2,3>*>(&rCellPopulation);


        // Calculate midpoint
        c_vector<double,3> centroid = zero_vector<double>(3);
        for (AbstractCellPopulation<2,3>::Iterator cell_iter = rCellPopulation.Begin();
                     cell_iter != rCellPopulation.End();
                     ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            centroid += rCellPopulation.GetNode(node_index)->rGetLocation();
        }
        centroid /= rCellPopulation.GetNumRealCells();

        for (AbstractCellPopulation<2,3>::Iterator cell_iter = rCellPopulation.Begin();
                     cell_iter != rCellPopulation.End();
                     ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<3>* p_node = rCellPopulation.GetNode(node_index);
            c_vector<double,3> cell_location = p_node->rGetLocation() - centroid;
            cell_location(2) = 0.0;
            c_vector<double, 3> force = zero_vector<double>(3);

            //Calculate cell normal (average of element normals)
            c_vector<double,3> normal = zero_vector<double>(3);

            std::set<unsigned>&  containing_elements = p_node->rGetContainingElementIndices();
            assert(containing_elements.size()>0);
            for (std::set<unsigned>::iterator iter = containing_elements.begin();
                    iter != containing_elements.end();
                    ++iter)
            {
                // Negative as normals point inwards for these surface meshes
                normal += - p_cell_population->rGetMesh().GetElement(*iter)->CalculateNormal();
            }
            normal /= norm_2(normal);

                

            double cell_area = rCellPopulation.GetVolumeOfCell(*cell_iter);
            force = mStrength * cell_area * normal; // cell_location / norm_2(cell_location);
            cell_iter->GetCellData()->SetItem("area", cell_area);
            
            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(force);

            cell_iter->GetCellData()->SetItem("applied_force_x", force[0]);
            cell_iter->GetCellData()->SetItem("applied_force_y", force[1]);
            cell_iter->GetCellData()->SetItem("applied_force_z", force[2]);
            cell_iter->GetCellData()->SetItem("applied_force_mag", norm_2(force));

            cell_iter->GetCellData()->SetItem("norm_x", normal[0]);
            cell_iter->GetCellData()->SetItem("norm_y", normal[1]);
            cell_iter->GetCellData()->SetItem("norm_z", normal[2]);

            cell_iter->GetCellData()->SetItem("radius", norm_2(cell_location));
        }
    }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<Strength>" << mStrength << "</Strength>\n";
        AbstractForce<2,3>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(RadialForce)

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(RadialForce)

class TestAnalyticComparisonLiteratePaper : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:

    void TestCylinderEquilateralImposedPressure() throw (Exception)
    {

        unsigned N_D[4] = {20,40,80,160};
        unsigned N_Z[6] = {15,30,60,120,240,480};

        for (unsigned N_D_index = 0; N_D_index < 4; N_D_index++)
        {
            for (unsigned N_Z_index = 0; N_Z_index < 3; N_Z_index++)
            {
                Honeycomb3DCylinderMeshGenerator generator(N_D[N_D_index], N_Z[N_D_index+N_Z_index], 1.5e-3, 12e-3);
                MutableMesh<2,3>* p_mesh = generator.GetMesh();
                p_mesh->Translate(-6.0e-3*unit_vector<double>(3,2));

                std::stringstream out;
                out << N_D[N_D_index] << "_" << N_Z[N_Z_index+N_D_index];
                std::string mesh_size = out.str();
                std::string output_directory = "CylinderValidation/Equilateral/" + mesh_size;

                // Create cells
                MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
                std::vector<CellPtr> cells;
                CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
                cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

                // Create a cell population
                MeshBasedCellPopulation<2,3> cell_population(*p_mesh, cells);
                cell_population.SetWriteVtkAsPoints(true);
                cell_population.SetOutputMeshInVtk(true);

                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellProliferativeTypesWriter>();
                cell_population.AddCellWriter<CellMutationStatesWriter>();

                cell_population.CalculateRestLengths();

                // Set up cell-based simulation
                OffLatticeSimulation<2,3> simulator(cell_population);
                simulator.SetOutputDirectory(output_directory);
                simulator.SetEndTime(M_TIME_FOR_SIMULATION); //50
                simulator.SetDt(0.002);
                simulator.SetSamplingTimestepMultiple(500);
                simulator.SetUpdateCellPopulationRule(false); // No remeshing.

                // Create a force law and pass it to the simulation
                boost::shared_ptr<GeneralisedLinearSpringForce<2,3> > p_linear_force(new GeneralisedLinearSpringForce<2,3>());
                p_linear_force->SetMeinekeSpringStiffness(50.0);

                simulator.AddForce(p_linear_force);

                // Create a force law to apply radial pressure force
                double pressure = 1.0666e4; // to match 80mmhg

                MAKE_PTR_ARGS(RadialForce, p_radial_force, (pressure));
                simulator.AddForce(p_radial_force);

                // Create a plane boundary to represent the inlet and pass them to the simulation
                boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition_1(new FixedRegionBoundaryCondition<2,3>(&cell_population,5.0e-3*unit_vector<double>(3,2),unit_vector<double>(3,2),10));
                simulator.AddCellPopulationBoundaryCondition(p_condition_1);

                boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition_2(new FixedRegionBoundaryCondition<2,3>(&cell_population,5.0e-3*unit_vector<double>(3,2),-unit_vector<double>(3,2),10));
                simulator.AddCellPopulationBoundaryCondition(p_condition_2);

                simulator.Solve();

                // To reset before looping: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
            }
        }
    }

void TestCylinderImposedPressureWithRandomMeshes() throw (Exception)
    {

        unsigned N[4] = {581,2103,8334,33331};
        
        for (unsigned N_index = 0; N_index < 4; N_index++)
        {
            std::stringstream out;
            out << N[N_index];
            std::string mesh_size = out.str();
            std::string mesh_file = "projects/EMBC2018/test/data/cyl_" + mesh_size + "_nodes.vtu";
            std::string output_directory = "CylinderValidation/Random/" + mesh_size;

            // This data file is in mm
            VtkMeshReader<2,3> mesh_reader(mesh_file);
            MutableMesh<2,3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            double scaling = 1e-3;  // so distances are in m

            mesh.Scale(scaling,scaling,scaling);

            // Create cells
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

            // Create a cell population
            MeshBasedCellPopulation<2,3> cell_population(mesh, cells);
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.SetOutputMeshInVtk(true);

            cell_population.AddCellWriter<CellIdWriter>();
            cell_population.AddCellWriter<CellProliferativeTypesWriter>();
            cell_population.AddCellWriter<CellMutationStatesWriter>();

            cell_population.CalculateRestLengths();

            // Set up cell-based simulation
            OffLatticeSimulation<2,3> simulator(cell_population);
            simulator.SetOutputDirectory(output_directory);
            simulator.SetEndTime(M_TIME_FOR_SIMULATION);
            simulator.SetDt(0.002);
            simulator.SetSamplingTimestepMultiple(500);
            simulator.SetUpdateCellPopulationRule(false); // No remeshing.

            // Create a force law and pass it to the simulation
            boost::shared_ptr<GeneralisedLinearSpringForce<2,3> > p_linear_force(new GeneralisedLinearSpringForce<2,3>());
            p_linear_force->SetMeinekeSpringStiffness(50.0);
            simulator.AddForce(p_linear_force);

            // Create a force law to apply radial pressure force
            double pressure = 1.0666e4; // to match 80mmhg

            MAKE_PTR_ARGS(RadialForce, p_radial_force, (pressure));
            simulator.AddForce(p_radial_force);

            // Create a plane boundary to represent the inlet and pass them to the simulation
            boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition_1(new FixedRegionBoundaryCondition<2,3>(&cell_population,scaling*5.0*unit_vector<double>(3,2),unit_vector<double>(3,2),10));
            simulator.AddCellPopulationBoundaryCondition(p_condition_1);

            boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition_2(new FixedRegionBoundaryCondition<2,3>(&cell_population,-scaling*5.0*unit_vector<double>(3,2),-unit_vector<double>(3,2),10));
            simulator.AddCellPopulationBoundaryCondition(p_condition_2);

            simulator.Solve();

            // To reset before looping: this is usually done by the SetUp and TearDown methods
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
        }
    }

};

#endif /*TESTCYLINDERVALIDATION_HPP_*/

