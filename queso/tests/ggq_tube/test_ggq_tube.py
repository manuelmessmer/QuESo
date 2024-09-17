# Project imports
from QuESo_PythonApplication.PyQuESo import PyQuESo

import unittest
import json

class TestGGQTube(unittest.TestCase):
    #p=2
    def test_1(self):
        self.RunTest("queso/tests/ggq_tube/QuESoSettings1.json", "queso/tests/ggq_tube/result_ips_gauss.json")

    def test_2(self):
        self.RunTest("queso/tests/ggq_tube/QuESoSettings2.json", "queso/tests/ggq_tube/result_ips_gauss_reduced1.json")

    def test_3(self):
        self.RunTest("queso/tests/ggq_tube/QuESoSettings3.json", "queso/tests/ggq_tube/result_ips_gauss_reduced2.json")

    def test_4(self):
        self.RunTest("queso/tests/ggq_tube/QuESoSettings4.json", "queso/tests/ggq_tube/result_ips_reduced_exact.json")

    def test_5(self):
        self.RunTest("queso/tests/ggq_tube/QuESoSettings5.json", "queso/tests/ggq_tube/result_ips_reduced_order1.json")

    def test_6(self):
        self.RunTest("queso/tests/ggq_tube/QuESoSettings6.json", "queso/tests/ggq_tube/result_ips_reduced_order2.json")

    def RunTest(self,filename, filename_result):
        self.pyqueso = PyQuESo(filename)
        self.pyqueso.Run()

        ips = {}
        num_el_inside = 0
        for i, element in enumerate(self.pyqueso.GetElements()):
            if not element.IsTrimmed():
                num_el_inside += 1
                tmp_list = []
                for point in element.GetIntegrationPoints():
                    tmp_list.append( [point.X(), point.Y(), point.Z(), point.Weight()] )
                ips[element.ID()] = tmp_list

        # with open("test.json", 'w') as file:
        #     json.dump(ips, file )

        ref_file = open(filename_result, "r")
        ref_ips = json.load(ref_file)
        ref_file.close()

        for ref_id, ref_el in ref_ips.items():
            el = ips[int(ref_id)]
            for ref_point, point in zip(ref_el, el):
                self.assertLess( abs(ref_point[0]-point[0])/ref_point[0], 1e-11 )
                self.assertLess( abs(ref_point[1]-point[1])/ref_point[1], 1e-11 )
                self.assertLess( abs(ref_point[2]-point[2])/ref_point[2], 1e-11 )
                self.assertLess( abs(ref_point[3]-point[3])/ref_point[3], 1e-11 )


if __name__ == "__main__":
    unittest.main()
