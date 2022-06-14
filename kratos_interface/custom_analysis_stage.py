# Import Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.gid_output_process import GiDOutputProcess

class CustomAnalysisStage(StructuralMechanicsAnalysis):
    def __init__(self, model, project_parameters, filename, volume_ref, lower_point, upper_point, boundary_conditions, integration_points_embedder, material_props, postprocess_flag, gid_output_dest):
        super().__init__(model, project_parameters)
        self.filename = filename
        self.volume_ref = volume_ref
        self.lower_point = lower_point
        self.upper_point = upper_point
        self.boundary_conditions = boundary_conditions
        self.gid_output_dest = gid_output_dest
        self.material_props = material_props
        tet_model_part = self.model.CreateModelPart("IgaModelPart")
        tet_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        tet_model_part.AddNodalSolutionStepVariable(KM.REACTION)
        self.integration_points_embedder = integration_points_embedder

        tet_model_part = self.model.CreateModelPart("TetModelPart")
        tet_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        tet_model_part.AddNodalSolutionStepVariable(KM.REACTION)

        self.Initialized = False
        self.time_ = []
        self.displacement = []
        self.analytical = []
        # Read model part
        if postprocess_flag:
            model_part_io_structure = KM.ModelPartIO(filename)
            model_part_io_structure.ReadModelPart(tet_model_part)

            KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, tet_model_part)
            KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, tet_model_part)
            KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, tet_model_part)
        # Read model part
        #model_part_io_structure = KM.ModelPartIO(self.filename)
        #model_part_io_structure.ReadModelPart(tet_model_part)

    def _ModelersSetupModelPart(self):
        model_part = self.model.GetModelPart('NurbsMesh')
        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.REACTION)

        model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE,3)

        # Create constitutive law and properties for BSpline elements
        cl = StructuralMechanicsApplication.LinearElastic3DLaw()
        volume_properties = model_part.GetProperties()[0]
        volume_properties.SetValue(KM.YOUNG_MODULUS, self.material_props[0])
        volume_properties.SetValue(KM.POISSON_RATIO, self.material_props[1])
        volume_properties.SetValue(KM.DENSITY, self.material_props[2])
        volume_properties.SetValue(KM.CONSTITUTIVE_LAW, cl)
        volume_properties.SetValue(KM.COMPUTE_LUMPED_MASS_MATRIX, False)

        # Convert the geometry model or import analysis suitable models.
        for modeler in self._GetListOfModelers():
            print(modeler)
            if self.echo_level > 1:
                KM.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup ModelPart started.")
            modeler.SetupModelPart()
            if self.echo_level > 1:
                KM.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup ModelPart finished.")

        print(model_part)
        nurbs_volume = model_part.GetGeometry("NurbsVolume")
        integration_points = []
        # Loop over all given integration points
        for index, point in enumerate(self.integration_points_embedder):
            integration_points.append([point.GetX(), point.GetY(), point.GetZ(), point.GetWeight()])

        # Create quadrature_point_geometries
        # These are basically just integration points.
        quadrature_point_geometries = KM.GeometriesVector()
        nurbs_volume.CreateQuadraturePointGeometries(quadrature_point_geometries, 2,
           integration_points)
        print("Number quadrature point geometries: ", len(quadrature_point_geometries))

        # Assing element formulation and constitutive law to each integration point.
        for i in range(0, len(quadrature_point_geometries)):
            el = model_part.CreateNewElement('SmallDisplacementElement3D8N', i+1, quadrature_point_geometries[i], volume_properties)

        print("Number of Elements/Integration Points (In BSpline Volume): ", model_part.NumberOfElements())

        # Add Dofs
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, model_part)

        return super()._ModelersSetupModelPart()

    def InitializeSolutionStep(self):
        if not self.Initialized:
            model_part = self.model.GetModelPart("NurbsMesh")
            for bc in self.boundary_conditions:
                bc.apply(model_part)
            self.Initialized = True

        return super().InitializeSolutionStep()

    def Timoshenko(self, t):
        import numpy as np
        return np.cos(t*1.522234e+3)*2.13715e-4+np.cos(t*8.843256e+1)*5.09311e-2+np.cos(t*5.542205e+2)*1.41972e-3-5.256458e-2

    def FinalizeSolutionStep(self):

        # tet_model_part = self.model.GetModelPart("IgaModelPart")
        # KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, tet_model_part)
        # KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, tet_model_part)
        # KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, tet_model_part)

        model_part = self.model.GetModelPart("NurbsMesh")
        nurbs_volume = model_part.GetGeometry("NurbsVolume")

        param = KM.Vector(3)
        param[0] = (0.0 - self.lower_point[0])/ abs(self.lower_point[0] - self.upper_point[0])
        param[1] = (0.0 - self.lower_point[1])/ abs(self.lower_point[1] - self.upper_point[1])
        param[2] = (10.0 - self.lower_point[2])/ abs(self.lower_point[2] - self.upper_point[2])

        self.time_.append(self.time)
        self.displacement.append( nurbs_volume.GlobalCoordinates(param)[2]-10 )

        # if( self.time < 3.000001 and self.time > 2.999999):
        #     local_disp = []
        #     import numpy as np
        #     x_values = np.arange(0,10.001,0.1)
        #     for x in x_values:
        #         param[2] = (x - self.lower_point[2])/ abs(self.lower_point[2] - self.upper_point[2])
        #         local_disp.append( (nurbs_volume.GlobalCoordinates(param)[2]-x) )

        #     import scipy
        #     from scipy.io import savemat
        #     mdic = {"x_coord": x_values, "y_coord": local_disp}
        #     savemat("snapshot.mat", mdic)
        #self.analytical.append( self.Timoshenko(self.time) )
        #self.displacement.append( -nurbs_volume.GlobalCoordinates(param)[1] )

        # for node in tet_model_part.Nodes:
        #     param = KM.Vector(3)
        #     param[0] = (node.X0 - self.lower_point[0])/abs(self.upper_point[0]-self.lower_point[0])
        #     param[1] = (node.Y0 - self.lower_point[1])/abs(self.upper_point[1]-self.lower_point[1])
        #     param[2] = (node.Z0 - self.lower_point[2])/abs(self.upper_point[2]-self.lower_point[2])
        #     # Check if params are inside bounds
        #     # if( param[0] < 0.0 or param[0] > 1.0):
        #     #     raise Exception("Analysis :: Postprocessing :: Parameter u out of bounds!")
        #     # if( param[1] < 0.0 or param[1] > 1.0):
        #     #     raise Exception("Analysis :: Postprocessing :: Parameter v out of bounds!")
        #     # if( param[2] < 0.0 or param[2] > 1.0):
        #     #     raise Exception("Analysis :: Postprocessing :: Parameter w out of bounds!")
        #     # Compute dispalcement with repsect to reference solution.
        #     coord_ref = self.volume_ref.GlobalCoordinates(param)
        #     coord = nurbs_volume.GlobalCoordinates(param)
        #     disp = coord - coord_ref
        #     # step = self._GetSolver().GetComputingModelPart().ProcessInfo[KM.STEP]
        #     # time = self._GetSolver().GetComputingModelPart().ProcessInfo[KM.TIME]

        #     # tet_model_part.ProcessInfo[KM.STEP] = step
        #     # tet_model_part.CloneTimeStep(time)

        #     node.SetSolutionStepValue(KM.DISPLACEMENT, 1, disp)



        return super().FinalizeSolutionStep()

    def Finalize(self):
        /import matplotlib.pyplot as plt


        # import scipy
        # from scipy.io import savemat
        # mdic = {"x_sim": self.time_, "y_sim": self.displacement}
        # savemat("matlab_mdisp.mat", mdic)
        return super().Finalize()

    def postprocess_mdpa(self, filename, nurbs_volume, volume_ref):
        print("postprocess_mdpa: ")
        # GiD Postprocessing
        ########################################################

        tet_model_part = self.model.GetModelPart("TetModelPart")
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, tet_model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, tet_model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, tet_model_part)
        tet_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        tet_model_part.AddNodalSolutionStepVariable(KM.REACTION)

        # #Read model part
        # model_part_io_structure = KM.ModelPartIO(filename)
        # model_part_io_structure.ReadModelPart(tet_model_part)

        KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, tet_model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, tet_model_part)
        KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, tet_model_part)

        for node in tet_model_part.Nodes:
            param = KM.Vector(3)
            param[0] = (node.X0 - self.lower_point[0])/abs(self.upper_point[0]-self.lower_point[0])
            param[1] = (node.Y0 - self.lower_point[1])/abs(self.upper_point[1]-self.lower_point[1])
            param[2] = (node.Z0 - self.lower_point[2])/abs(self.upper_point[2]-self.lower_point[2])
            # Check if params are inside bounds
            # if( param[0] < 0.0 or param[0] > 1.0):
            #     raise Exception("Analysis :: Postprocessing :: Parameter u out of bounds!")
            # if( param[1] < 0.0 or param[1] > 1.0):
            #     raise Exception("Analysis :: Postprocessing :: Parameter v out of bounds!")
            # if( param[2] < 0.0 or param[2] > 1.0):
            #     raise Exception("Analysis :: Postprocessing :: Parameter w out of bounds!")
            # Compute dispalcement with repsect to reference solution.
            coord_ref = self.volume_ref.GlobalCoordinates(param)
            coord = nurbs_volume.GlobalCoordinates(param)
            disp = coord - coord_ref
            # step = self._GetSolver().GetComputingModelPart().ProcessInfo[KM.STEP]
            # time = self._GetSolver().GetComputingModelPart().ProcessInfo[KM.TIME]

            # tet_model_part.ProcessInfo[KM.STEP] = step
            # tet_model_part.CloneTimeStep(time)

            node.SetSolutionStepValue(KM.DISPLACEMENT, 1, disp)

        # KM.VariableUtils().AddDof(KM.DISPLACEMENT_X, KM.REACTION_X, tet_model_part)
        # KM.VariableUtils().AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y, tet_model_part)
        # KM.VariableUtils().AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z, tet_model_part)

        #         node.SetSolutionStepValue(KM.DISPLACEMENT,disp)
        # import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle
        # from KratosMultiphysics.ParticleMechanicsApplication.particle_gid_output_process import ParticleGiDOutputProcess
        output_file = self.gid_output_dest + "gid_output_"
        gid_output =  GiDOutputProcess(tet_model_part,
                                output_file,
                                KM.Parameters("""
                                    {
                                        "result_file_configuration": {
                                            "gidpost_flags": {
                                                "GiDPostMode": "GiD_PostAscii",
                                                "WriteDeformedMeshFlag": "WriteUndeformed",
                                                "MultiFileFlag": "SingleFile"
                                            },
                                            "file_label": "step",
                                            "output_control_type": "step",
                                            "output_interval": 1.0,
                                            "body_output": true,
                                            "node_output": false,
                                            "skin_output": false,
                                            "plane_output": [],
                                            "nodal_results": ["DISPLACEMENT"],
                                            "nodal_nonhistorical_results": [],
                                            "nodal_flags_results": [],
                                            "gauss_point_results": [],
                                            "additional_list_files": []
                                        }
                                    }
                                    """)
                            )

        gid_output.ExecuteInitialize()
        gid_output.ExecuteBeforeSolutionLoop()
        gid_output.ExecuteInitializeSolutionStep()
        gid_output.PrintOutput()
        gid_output.ExecuteFinalizeSolutionStep()
        gid_output.ExecuteFinalize()