import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

"""
This is an example of dummy_solver.py.
To use test functions for the testing of a single solver, a file with this name should be included in the working directory.
The functions that need to be defined depend on the variables for the interface.
The names of these functions are fixed as calculate_<variable>(x,y,z,n).
The functions receive the x, y and z-coordinate of the nodes in undeformed state and the current time step (n).
They have to return a list or numpy array of 1 or 3 elements for a scalar or vector, respectively.
Several types of test can be grouped into this dummy_solver.py file by creating additional classes.
The name of the class to be used should be specified in the .json file containing the settings for the case.
"""


class SimpleTest:
    def calculate_displacement(self, x, y, z, n):
        """ Specify the displacement of a point on the interface relative to its original position"""
        if 0.01 < x < 0.024:
            disp = [0, 0, 0.001]
        else:
            disp = [0, 0, 0]
        return disp

    def calculate_pressure(self, x, y, z, n):
        """ Specify the pressure on the surface"""
        if x > 0.01:
            pres = [1000]
        else:
            pres = [0]
        return pres

    def calculate_traction(self, x, y, z, n):
        """ Specify the traction on the surface"""
        if x > 0.01:
            trac = [0, 0, 100]
        else:
            trac = [0, 0, 0]
        return trac


class TransientTest:
    def calculate_displacement(self, x, y, z, n):
        """ Specify the displacement of a point on the interface relative to its original position"""
        if n < 5:
            if 0.01 < x < 0.024:
                disp = [0, 0, 0.001]
            else:
                disp = [0, 0, 0]
        else:
            if 0.01 < x < 0.024:
                disp = [0, 0, -0.001]
            else:
                disp = [0, 0, 0]
        return disp

    def calculate_pressure(self, x, y, z, n):
        """ Specify the pressure on the surface"""
        if n < 5:
            if x > 0.01:
                pres = [1000]
            else:
                pres = [0]
        else:
            if x > 0.01:
                pres = [-1000]
            else:
                pres = [0]
        return pres

    def calculate_traction(self, x, y, z, n):
        """ Specify the traction on the surface"""
        if n < 5:
            if x > 0.01:
                trac = [0, 0, 100]
            else:
                trac = [0, 0, 0]
        else:
            if x > 0.01:
                trac = [0, 100, 0]
            else:
                trac = [0, 0, 0]
        return trac


class PythonSolverTest:
    def calculate_displacement(self, x, y, z, n):
        """ Specify the displacement of a point on the interface relative to its original position"""
        if n < 5:
            if 0.01 < z < 0.024:
                disp = [0, 0.0002, 0]
            else:
                disp = [0, 0, 0]
        else:
            if 0.01 < z < 0.024:
                disp = [0, -0.0001, 0]
            else:
                disp = [0, 0, 0]
        return disp

    def calculate_pressure(self, x, y, z, n):
        """ Specify the pressure on the surface"""
        if n < 5:
            if z > 0.01:
                pres = [1000]
            else:
                pres = [0]
        else:
            if z > 0.01:
                pres = [-1000]
            else:
                pres = [0]
        return pres

    def calculate_traction(self, x, y, z, n):
        """ Specify the traction on the surface"""
        if n < 5:
            if z > 0.01:
                trac = [0, 0, 100]
            else:
                trac = [0, 0, 0]
        else:
            if z > 0.01:
                trac = [0, 100, 0]
            else:
                trac = [0, 0, 0]
        return trac


class InterpolatedData:
    def __init__(self):
        abscissa = np.array([-0.025, -0.02433333, -0.02366667, -0.023, -0.02233333,
                             -0.02166667, -0.021, -0.02, -0.01666667, -0.01333333,
                             -0.01, -0.00666667, -0.00333333, 0., 0.00333333,
                             0.00666667, 0.01, 0.01333333, 0.01666667, 0.02,
                             0.02066667, 0.02133333, 0.022, 0.02266667, 0.02333333,
                             0.024, 0.02466667])
        y_displacement = np.array([1.16621303e-10, 1.95088490e-04, 3.60966621e-04, 4.83640889e-04,
                                   5.66206281e-04, 6.17896253e-04, 6.48283040e-04, 6.70057566e-04,
                                   6.78144584e-04, 6.75857193e-04, 6.74728085e-04, 6.73872500e-04,
                                   6.73091675e-04, 6.72355520e-04, 6.71655600e-04, 6.70986983e-04,
                                   6.70345131e-04, 6.69822764e-04, 6.70441553e-04, 6.63535153e-04,
                                   6.52518454e-04, 6.31071264e-04, 5.92546400e-04, 5.28036801e-04,
                                   4.27512507e-04, 2.83584753e-04, 9.96992219e-05])
        self.interpolant = interpolate.splrep(abscissa, y_displacement, s=0)
        plt.plot(abscissa, y_displacement)

    def calculate_displacement(self, x, y, z, n):

        #""" Specify the displacement of a point on the interface relative to its original position"""
        if -0.025 <= x <= 0.025:
            disp = [0, float(interpolate.splev(x, self.interpolant, der=0)), 0]
        else:
            disp = [0, 0, 0]
        return disp

class boundaryData:
    def __init__(self):
        points = np.array([-0.00296500005411233, -0.00289499991413968, -0.00282500033512334, -0.00275500024178281,
                           -0.00268500022903399, -0.00261500030262096,
                           -0.00254500032142289, -0.00247500030084457, -0.00240500022992312, -0.00233500023945582,
                           -0.00226500027045965, -0.00219500028623788,
                           -0.00212500027873324, -0.00205500028147257, -0.00198500029107292, -0.00191500029001559,
                           -0.00184500029731382, -0.00177500028688709,
                           -0.00170500027288193, -0.00163500024983635, -0.00156500023110855, -0.00149500026206993,
                           -0.00142500025728856, -0.00135500028595614,
                           -0.00128500027552768, -0.00121500029031525, -0.00114500029880467, -0.00107500027977865,
                           -0.00100500027711092, -0.000935000314155058,
                           -0.000865000342132312, -0.000795000247870493, -0.000725000271622883, -0.000655000085681139,
                           -0.000585000305934551, -0.000515000173635347,
                           -0.000444997288207446, -0.000375002356527749, -0.000305002962168476, -0.000235002959938602,
                           -0.000165002980697297, -9.50030304288466e-05,
                           -2.50029796287272e-05, 4.49970272374523e-05, 0.000114997061919652, 0.000184997171976562,
                           0.000254997020037823, 0.000324997012247883,
                           0.000395000352129606, 0.000464997503757271])

        pressure = np.array([
            2363098.62428648, 7283648.24537446, 12498383.7888241, 17961389.1754059, 23713543.9629163, 29778639.4177025,
            36177428.2057563,
            42937440.9688923, 50090048.0132484, 57674651.9161937, 65729348.3645303, 74296340.2261822, 83424459.8155285,
            93170917.9794775,
            103599443.612899, 114782032.499944, 126802795.014243, 139756842.158668, 153756750.428036, 168933593.9815,
            185442998.100372,
            203469847.703922, 223222989.548409, 244962976.879889, 268993107.885952, 295696009.664748, 325533669.162073,
            359080796.299739,
            397075217.42787, 440460606.093902, 490431750.202189, 548546432.955527, 617104447.319751, 699126107.457369,
            799291148.44732,
            923350791.466108, 1080156299.8478, 1304505005.3195, 1606624529.10661, 1996515658.29936, 2601125208.37708,
            3345115953.42972,
            5244651206.39622, 4687269103.86258, 4762840107.35952, 4034559510.34207, 3305402301.49215, 2367644874.79225,
            1428726372.25248,
            478799253.5599])
        self.interpolant_pressure = interpolate.interp1d(points, pressure, fill_value="extrapolate")
        plt.plot(points,pressure)

        traction = np.array([[7845797.41728537, - 45110.4473901315, 1.05937730117891e-10],
                                 [8001786.2159395, - 33321.0638375707, 8.47634421396717e-10],
                                 [8198483.80279624, - 36350.0723190803, 1.52203421519673e-09],
                                 [8397428.01259482, - 39013.3006409989, 1.55541901087833e-09],
                                 [8612205.69946456, - 40025.2908848634, 1.80001993378652e-10],
                                 [8837444.75754526, - 41414.796338382, 9.54285646943498e-10],
                                 [9073218.68853115, - 42773.3802060612, - 2.15979714467827e-09],
                                 [9321996.09203247, - 43623.802688657, - 1.14438388506531e-09],
                                 [9586434.45770407, - 44436.8516397956, - 8.49398988461123e-10],
                                 [9867390.53771502, - 45600.5689826222, - 3.17157930573979e-11],
                                 [10164364.3188675, - 47194.1380114693, 1.07952453071507e-09],
                                 [10478895.6920955, - 48745.1551251216, - 4.59130403481607e-09],
                                 [10813536.8835585, - 50211.3794694182, - 4.39778155898744e-09],
                                 [11170485.5031434, - 51777.0285273675, - 1.78218007136222e-10],
                                 [11551515.1926222, - 53483.4000966924, - 3.64837848229356e-09],
                                 [11959301.6880953, - 55273.7759536723, 1.95554000601409e-09],
                                 [12396780.6203422, - 57129.8327015677, - 1.18163093376948e-09],
                                 [12867501.5235592, - 58985.7070310552, - 5.04197353864074e-09],
                                 [13375975.1018488, - 60819.115375995, - 3.90042509208466e-09],
                                 [13927302.5720987, - 62658.3879483826, 5.16606057245274e-09],
                                 [14527449.8080229, - 64859.1665550307, - 8.10762803260661e-11],
                                 [15181071.1242984, - 67635.0681756715, 5.18041794134327e-09],
                                 [15895290.2883942, - 70524.3616528996, - 5.86621200860734e-09],
                                 [16679092.8115055, - 73715.2299749723, - 1.19246286214699e-08],
                                 [17543470.8998415, - 76837.7264959746, - 3.00541155127441e-09],
                                 [18502241.2721205, - 80438.7854740059, 2.25634028500227e-10],
                                 [19570287.2815257, - 84234.1616333273, 5.51616150276804e-09],
                                 [20769878.2624432, - 87716.9619752825, 2.94069221547837e-09],
                                 [22128176.6927192, - 91862.9384469539, 6.11058559430309e-11],
                                 [23673070.6454734, - 97598.5896368533, 4.74447530911987e-09],
                                 [25442730.2188137, - 103402.192228648, - 9.84584889698139e-09],
                                 [27507920.9546045, - 105992.347098972, - 4.15338974565324e-10],
                                 [29942024.5077651, - 112350.531748138, - 1.0421826268095e-08],
                                 [32876316.9735011, - 106651.109816807, 1.83064317838841e-08],
                                 [36427907.8368433, - 137587.005851713, 3.33119408537761e-10],
                                 [40787841.0715227, - 178000.571855228, 9.30531448934348e-09],
                                 [47339678.0659371, 346508.526343755, 4.75402864955976e-09],
                                 [55716795.0386784, - 533822.062750933, 1.42005563969614e-08],
                                 [63261017.0455455, - 350073.917304143, - 3.41296819172703e-09],
                                 [75096592.4818617, - 292140.858013722, - 2.37987470770375e-09],
                                 [84652006.2405111, - 83135.1810068885, 8.0189969355175e-09],
                                 [114550331.1408, 83387.3338945355, 1.56662093175585e-08],
                                 [90493101.2013239, 107664.740234037, 2.03796597497221e-08],
                                 [65711046.4252014, - 229887.373760638, 7.18180316446074e-09],
                                 [56535339.3540119, - 227337.860533094, 2.41251481352112e-09],
                                 [36526893.0814791, - 622132.094404016, 1.08420938550205e-08],
                                 [26036308.2675188, - 361242.599498109, - 9.43012577180967e-09],
                                 [15176530.1673735, 3365410.88619825, - 5.99094884886286e-09],
                                 [12407094.0802307, -7012396.84468801, -5.11367251201859e-09],
                                 [9584855.25094096, 8256248.01610988, - 8.55696100428801e-12]])


        self.interpolant_traction_0 = interpolate.interp1d(points, traction[:,0])
        self.interpolant_traction_1 = interpolate.interp1d(points, traction[:,1])
        self.interpolant_traction_2 = interpolate.interp1d(points, traction[:,2])

    def calculate_pressure(self, x, y, z, n):
        print(x)
        pressure = [self.interpolant_pressure(x)]
        print(pressure)
        return pressure

    def calculate_traction(self, x, y, z, n):
        if -0.003 <= x <= 0.0005:
            traction  = [self.interpolant_traction_0(x), self.interpolant_traction_1(x), self.interpolant_traction_0[2]]
        return traction
