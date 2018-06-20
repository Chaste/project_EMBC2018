#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

#include <cxxtest/TestSuite.h>
#include <cstdio>
#include <ctime>
#include <cmath>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "OffLatticeSimulation.hpp"
#include "AppliedForceModifier.hpp"

#include "CellIdWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"

#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "AppliedForce.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "SmartPointers.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VtkMeshWriter.hpp"
#include "CommandLineArguments.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"
#include "XmlTools.hpp"
using namespace xsd::cxx::tree;


class TestSetupFlowInPipe : public AbstractCellBasedTestSuite
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

    void ReadHemeLBv3Vector(c_vector<double,3>& vector, xercesc::DOMElement* xmlElement)
    {
        std::stringstream raw_value(X2C(xmlElement->getAttribute(X("value"))));
        char left_par, comma1, comma2, right_par;
        double x, y, z;
        raw_value >> left_par; assert(left_par == '(');
        raw_value >> x;
        raw_value >> comma1; assert(comma1 == ',');
        raw_value >> y;
        raw_value >> comma2; assert(comma2 == ',');
        raw_value >> z;
        raw_value >> right_par; assert(right_par == ')');

        vector[0] = x;
        vector[1] = y;
        vector[2] = z;
    }

    void ReadIoletPlanes(std::string hemelbConfigFile,
                         std::vector<c_vector<double,3> >& rBoundaryPlanePoints,
                         std::vector<c_vector<double,3> >& rBoundaryPlaneNormals,
                         std::vector<double>& rBoundaryPlaneRadii)
    {
        PRINT_VARIABLE(hemelbConfigFile);
        bool hemelb_config_v3;
        std::vector<xercesc::DOMElement*> iolets;
        try
        {
        	properties<char> props;
            xsd::cxx::xml::dom::auto_ptr<xercesc::DOMDocument> p_doc = XmlTools::XmlTools::ReadXmlFile(hemelbConfigFile, props, false);
            xercesc::DOMElement* p_root_elt = p_doc->getDocumentElement();
            hemelb_config_v3 = (atoi(X2C(p_root_elt->getAttribute(X("version"))).c_str()) == 3);
            iolets = XmlTools::FindElements(p_root_elt, "inlets/inlet");
            std::vector<xercesc::DOMElement*> aux = XmlTools::FindElements(p_root_elt, "outlets/outlet");
            iolets.insert(iolets.end(), aux.begin(), aux.end());
        }
        catch (const exception<char>& e)
        {
            std::cerr << e << std::endl;
            EXCEPTION("XML parsing error in configuration file: " + hemelbConfigFile);
        }

        if (!hemelb_config_v3)
        {
            EXCEPTION("Only HemeLB XML version 3 is supported");
        }

        for (unsigned iolet_id=0; iolet_id<iolets.size(); iolet_id++)
        {
            c_vector<double,3> vector;

            std::vector<xercesc::DOMElement*> points = XmlTools::FindElements(iolets[iolet_id], "position");
            assert(points.size()==1);
			ReadHemeLBv3Vector(vector, points[0]);
            rBoundaryPlanePoints.push_back(vector);

            std::vector<xercesc::DOMElement*> normals = XmlTools::FindElements(iolets[iolet_id], "normal");
            assert(normals.size()==1);
			ReadHemeLBv3Vector(vector, normals[0]);
            rBoundaryPlaneNormals.push_back(vector);

            std::vector<xercesc::DOMElement*> radi = XmlTools::FindElements(iolets[iolet_id], "radius");
            if(radi.size()==1)
            {
            	double radius = atof(X2C(radi[0]->getAttribute(X("radius"))).c_str());
            	rBoundaryPlaneRadii.push_back(radius);
            }
            else
            {
            	rBoundaryPlaneRadii.push_back(40e-6); // in m
            }
        }
    }

public:

    void TestSetupRunAndSavePipe() throw (Exception)
    {
        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-mesh"));
        std::string mesh_file = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-mesh");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-mesh_scale"));
        double mesh_scale = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-mesh_scale").c_str());

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-xml"));
        std::string hemelb_config_file = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-xml");

        double spring_constant = 1.0;
        double membrane_constant = 1e-12;
        double drag_coeficient = 1.0;

        // Read inlet position from the Hemelb input XML file
        std::vector<c_vector<double,3> > boundary_plane_points;
        std::vector<c_vector<double,3> > boundary_plane_normals;
        std::vector<double> boundary_plane_radii;
        ReadIoletPlanes(hemelb_config_file, boundary_plane_points, boundary_plane_normals, boundary_plane_radii);

        VtkMeshReader<2,3> mesh_reader(mesh_file);
        MutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Scale(mesh_scale,mesh_scale,mesh_scale); // so distances are in m

        // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        MeshBasedCellPopulation<2,3> cell_population(mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        cell_population.CalculateRestLengths();
        cell_population.SetDampingConstantNormal(drag_coeficient);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2,3> simulator(cell_population);
        std::string output_dir = "./";
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetUpdateCellPopulationRule(false);

        // Create an Applied Force modifier to couple to Flow
        boost::shared_ptr<AppliedForceModifier<2,3> > p_force_modifier(new AppliedForceModifier<2,3>());
        p_force_modifier->SetResetTractionsOnCells(false,"");
        simulator.AddSimulationModifier(p_force_modifier);

        // Elastic force
        boost::shared_ptr<GeneralisedLinearSpringForce<2,3> > p_force ( new GeneralisedLinearSpringForce<2,3>() );
        p_force->SetMeinekeSpringStiffness(spring_constant);
        simulator.AddForce(p_force);

        // Create a plane boundary to represent the inlets/outlets and pass them to the simulation
        double inlet_offset = 0.0;
        for(unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
            boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population,boundary_plane_points[boundary_id]+inlet_offset*boundary_plane_normals[boundary_id],-boundary_plane_normals[boundary_id],boundary_plane_radii[boundary_id]));
            simulator.AddCellPopulationBoundaryCondition(p_condition);
        }

        // Save the set up simulation ready to be executed once flow results are available
        CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(&simulator);
    
        // Output the mesh in .vtu format for HemeLB setup tool to pick up (first converted to stl, though).
        VtkMeshWriter<2,3> mesh_writer(output_dir, "config", false);
        MutableMesh<2,3>* p_mesh = &(dynamic_cast<MeshBasedCellPopulation<2,3>*>(&(simulator.rGetCellPopulation()))->rGetMesh());
        p_mesh->Scale(1.0/mesh_scale,1.0/mesh_scale,1.0/mesh_scale); // so distances are back in original scale
        mesh_writer.WriteFilesUsingMesh(*p_mesh);
    }
};

#endif /*TESTRELAXATION_HPP_*/

