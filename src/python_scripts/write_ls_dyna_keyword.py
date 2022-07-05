import numpy as np

def truncate(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    n = n - len(i)
    return '.'.join([i, (d+'0'*n)[:n]])

class keyword_writer:
    def write_nurbs_volume(filename, volume):
        with open(filename + ".k","w") as file:
            file.write("*KEYWORD\n")
            file.write("*IGA_3D_NURBS_XYZ\n")
            file.write('{: <2}{: >8}{: >10}{: >10}{: >10}{: >10}{: >10}{: >10}\n'.format("$#", "nid", "nr", "ns", "nt", "pr", "ps", "pt") )
            file.write('{: >10}{: >10}{: >10}{: >10}{: >10}{: >10}{: >10}\n'.format(1, volume.NumberOfControlPointsU(),
                                                                                       volume.NumberOfControlPointsV(),
                                                                                       volume.NumberOfControlPointsW(),
                                                                                       volume.PolynomialDegreeU(),
                                                                                       volume.PolynomialDegreeV(),
                                                                                       volume.PolynomialDegreeW() ))
            file.write('{: <2}{: >8}{: >10}{: >10}\n'.format("$#", "unir", "unis", "unit") )
            file.write('{: >10}{: >10}{: >10}\n'.format(0,0,0) )
            #Knots in u
            file.write('{:<2}{: >18}{: >20}{: >20}{: >20}\n'.format("$#", "r1", "r2", "r3", "r4") )
            file.write('{: >20}'.format(volume.KnotsU()[0]) )
            for i, knot in enumerate(volume.KnotsU()):
                file.write('{: >20}'.format(knot) )
                if (i+2)%4==0:
                    file.write('\n')
            file.write('{: >20}'.format(volume.KnotsU()[len(volume.KnotsU())-1]) )
            file.write('\n')
            #Knots in v
            file.write('{:<2}{: >18}{: >20}{: >20}{: >20}\n'.format("$#", "s1", "s2", "s3", "s4") )
            file.write('{: >20}'.format(volume.KnotsV()[0]) )
            for i, knot in enumerate(volume.KnotsV()):
                file.write('{: >20}'.format(knot) )
                if (i+2)%4==0:
                    file.write('\n')
            file.write('{: >20}'.format(volume.KnotsV()[len(volume.KnotsV())-1]) )
            file.write('\n')
            # Knots in w
            file.write('{:<2}{: >18}{: >20}{: >20}{: >20}\n'.format("$#", "t1", "t2", "t3", "t4") )
            file.write('{: >20}'.format(volume.KnotsW()[0]) )
            for i, knot in enumerate(volume.KnotsW()):
                file.write('{: >20}'.format(knot) )
                if (i+2)%4==0:
                    file.write('\n')
            file.write('{: >20}'.format(volume.KnotsW()[len(volume.KnotsW())-1]) )
            file.write('\n')
            # Control points
            file.write('{:<2}{: >18}{: >20}{: >20}{: >20}\n'.format("$#", "x", "y", "z", "wgt") )
            for point in volume:
                x_out = point.X0
                y_out = point.Y0
                z_out = point.Z0
                if np.abs(x_out) < 1e-14:
                    x_out = 0.0
                if np.abs(y_out) < 1e-14:
                    y_out = 0.0
                if np.abs(z_out) < 1e-14:
                    z_out = 0.0
                file.write('{: >20}{: >20}{: >20}{: >20}\n'.format(truncate(x_out,16), truncate(y_out,16), truncate(z_out,16), 1.0) )
            file.write("*End")

    def write_points(filename, elements, ntot):
        to_write = {}
        total_number_points = 0
        total_number_elements = ntot
        for element in elements:
            array = []
            number_elements = 0
            if element.IsTrimmed():
                for point_trimmed in element.GetIntegrationPointsTrimmed():
                    if( point_trimmed.GetWeight() > 0.0 ):
                        number_elements -= 1
                        total_number_points += 1
                        xx = point_trimmed.GetX()
                        yy = point_trimmed.GetY()
                        zz = point_trimmed.GetZ()
                        ww = point_trimmed.GetWeight()
                        array.append( [xx, yy, zz, ww] )
            else:
                number_elements = len(element.GetIntegrationPointsInside())
                total_number_points += number_elements
                for point_inside in element.GetIntegrationPointsInside():
                    xx = point_inside.GetX()
                    yy = point_inside.GetY()
                    zz = point_inside.GetZ()
                    ww = point_inside.GetWeight()
                    array.append( [xx,yy,zz,ww] )

            to_write[element.ID()] = {'nip' : number_elements, 'ip' : array}

        with open(filename + ".k","w") as file:
            file.write("*KEYWORD\n")
            file.write("*IGA_USER_IPTS_SOLID\n")
            file.write('{: <2}{: >8}{: >10}{: >10}{: >10}{: >10}\n'.format("$#", "npid", "nel", "niptot", "pidis", "pidih") )
            file.write('{: >10}{: >10}{: >10}{: >10}{: >10}\n'.format(1, total_number_elements, total_number_points, '', 99) )
            file.write('{:}{:}\n'.format("$nip", " nip"*19 ) )
            for i in range(total_number_elements):
                if i+1 in to_write:
                    nip = to_write[i+1]["nip"]
                    file.write('{: >4}'.format(nip) )
                else:
                    file.write('{: >4}'.format(0) )
                if (i+1)%20 == 0:
                    file.write('\n')
            if( total_number_elements % 20 != 0):
                file.write('\n')
            file.write('{: <2}{: >18}{: >20}{: >20}{: >20}\n'.format("$#","ipr","ips", "ipt", "wgt") )
            for i in range(total_number_elements):
                if i+1 in to_write:
                    array = to_write[i+1]["ip"]
                    for ip in array:
                        file.write('{: >20.17f}{: >20.17f}{: >20.17f}{: >20.17f}\n'.format(ip[0], ip[1], ip[2], ip[3]) )
            file.write('*End')