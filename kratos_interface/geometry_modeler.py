import KratosMultiphysics
from KratosMultiphysics.modeler_factory import KratosModelerFactory
import KratosMultiphysics.IgaApplication as IgaApplication

class GeometryModeler:
    def run_modeler(current_model, modelers_list):
        factory = KratosModelerFactory()
        list_of_modelers = factory.ConstructListOfModelers(current_model, modelers_list)
        for modeler in list_of_modelers:
            modeler.SetupGeometryModel()

        for modeler in list_of_modelers:
            modeler.SetupModelPart()

    def create_geometry(model, lower_point, upper_point, number_of_knot_spans, polynomial_order):
        # Here, we create the background mesh
        # From point A to point B
        # number of knot spans equals number of elements
        modeler_settings = KratosMultiphysics.Parameters("""
        [{
            "modeler_name": "NurbsGeometryModeler",
            "Parameters": {
                "model_part_name" : "NurbsMesh",
                "lower_point": [-0.2, -0.2, -0.2],
                "upper_point": [0.2, 0.2, 0.4],
                "polynomial_order" : [2, 2, 2],
                "number_of_knot_spans" : [1,1,1]
            }
        }]
        """)

        # Override background grid parameters
        modeler_settings[0]["Parameters"]["lower_point"].SetVector( lower_point )
        modeler_settings[0]["Parameters"]["upper_point"].SetVector( upper_point )
        modeler_settings[0]["Parameters"]["polynomial_order"].SetVector( polynomial_order )
        modeler_settings[0]["Parameters"]["number_of_knot_spans"].SetVector( number_of_knot_spans )

        GeometryModeler.run_modeler(model, modeler_settings)

        return modeler_settings