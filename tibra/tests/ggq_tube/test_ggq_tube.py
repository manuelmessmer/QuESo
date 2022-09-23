# Project imports
from TIBRA_PythonApplication.PyTIBRA import PyTIBRA

import unittest
import json

class TestGGQTube(unittest.TestCase):
    #p=2
    def test_1(self):
        self.RunTest("tibra/tests/ggq_tube/TIBRAParameters1.json", "tibra/tests/ggq_tube/result_ips_gauss.json")

    def test_2(self):
        self.RunTest("tibra/tests/ggq_tube/TIBRAParameters2.json", "tibra/tests/ggq_tube/result_ips_gauss_reduced1.json")

    def test_3(self):
        self.RunTest("tibra/tests/ggq_tube/TIBRAParameters3.json", "tibra/tests/ggq_tube/result_ips_gauss_reduced2.json")

    def test_4(self):
        self.RunTest("tibra/tests/ggq_tube/TIBRAParameters4.json", "tibra/tests/ggq_tube/result_ips_reduced_exact.json")

    def test_5(self):
        self.RunTest("tibra/tests/ggq_tube/TIBRAParameters5.json", "tibra/tests/ggq_tube/result_ips_reduced_order1.json")

    def test_6(self):
        self.RunTest("tibra/tests/ggq_tube/TIBRAParameters6.json", "tibra/tests/ggq_tube/result_ips_reduced_order2.json")

    def RunTest(self,filename, filename_result):
        self.pytibra = PyTIBRA(filename)
        self.pytibra.Run()

        ips = {}
        num_el_inside = 0
        for i, element in enumerate(self.pytibra.GetElements()):
            if not element.IsTrimmed():
                num_el_inside += 1
                tmp_list = []
                for point in element.GetIntegrationPointsInside():
                    tmp_list.append( [point.GetX(), point.GetY(), point.GetZ(), point.GetWeight()] )
                ips[element.ID()] = tmp_list

        # with open("result_ips_gauss_reduced2.json", "w") as fp:
        #     json.dump(ips, fp)

        ref_file = open(filename_result, "r")
        ref_ips = json.load(ref_file)
        ref_file.close()

        for ref_id, ref_el in ref_ips.items():
            el = ips[int(ref_id)]
            for ref_point, point in zip(ref_el, el):
                self.assertLess( abs(ref_point[0]-point[0])/ref_point[0], 1e-14 )
                self.assertLess( abs(ref_point[1]-point[1])/ref_point[1], 1e-14 )
                self.assertLess( abs(ref_point[2]-point[2])/ref_point[2], 1e-14 )
                self.assertLess( abs(ref_point[3]-point[3])/ref_point[3], 1e-14 )


if __name__ == "__main__":
    unittest.main()
