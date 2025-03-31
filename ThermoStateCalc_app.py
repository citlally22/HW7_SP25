#region imports
import sys
from ThermoStateCalc import Ui__frm_StateCalculator # Import the UI generated by PyQt
from pyXSteam.XSteam import XSteam # Import pyXSteam for thermodynamic calculations
from PyQt5.QtWidgets import QWidget, QApplication # Import PyQt5 for GUI
from UnitConversion import UC # Import unit conversion utilities
from scipy.optimize import fsolve # Import fsolve for solving equations

# Class to handle saturation properties of a substance
class thermoSatProps:
    """
        A class to calculate saturation properties of a substance at a given pressure or temperature.

        Attributes:
            steamTable (XSteam): Instance of XSteam for thermodynamic calculations.
            pSat (float): Saturation pressure (bar).
            tSat (float): Saturation temperature (C).
            vf (float): Specific volume of saturated liquid (m^3/kg).
            vg (float): Specific volume of saturated vapor (m^3/kg).
            hf (float): Enthalpy of saturated liquid (kJ/kg).
            hg (float): Enthalpy of saturated vapor (kJ/kg).
            uf (float): Internal energy of saturated liquid (kJ/kg).
            ug (float): Internal energy of saturated vapor (kJ/kg).
            sf (float): Entropy of saturated liquid (kJ/kg·K).
            sg (float): Entropy of saturated vapor (kJ/kg·K).
            vgf (float): Difference in specific volume between vapor and liquid (m^3/kg).
            hgf (float): Difference in enthalpy between vapor and liquid (kJ/kg).
            sgf (float): Difference in entropy between vapor and liquid (kJ/kg·K).
            ugf (float): Difference in internal energy between vapor and liquid (kJ/kg).
        """
    def __init__(self, p=None, t=None, SI=True):
        """
        Initialize the thermoSatProps class.

        Args:
            p (float, optional): Pressure in bar. Defaults to None.
            t (float, optional): Temperature in Celsius. Defaults to None.
            SI (bool, optional): Use SI units if True, English units if False. Defaults to True.
        """
        # Initialize the steam table with the appropriate unit system
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS if SI else XSteam.UNIT_SYSTEM_FLS)
        # If pressure is provided, calculate saturation properties
        if p is not None:
            self.getSatProps(p, SI)
        # If temperature is provided, calculate saturation properties
        elif t is not None:
            self.getSatProps(self.steamTable.psat_t(t), SI)

    def getSatProps(self, p, SI=True):
        """
        Calculate saturation properties at the given pressure.

        Args:
            p (float): Pressure in bar.
            SI (bool, optional): Use SI units if True, English units if False. Defaults to True.
        """
        # Reinitialize steam table with the correct unit system
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS if SI else XSteam.UNIT_SYSTEM_FLS)
        # Store saturation pressure and calculate corresponding temperature
        self.pSat = p
        self.tSat = self.steamTable.tsat_p(p)
        # Calculate saturation properties for liquid and vapor
        self.vf = self.steamTable.vL_p(p)  # Specific volume of liquid
        self.vg = self.steamTable.vV_p(p)  # Specific volume of vapor
        self.hf = self.steamTable.hL_p(p)  # Enthalpy of liquid
        self.hg = self.steamTable.hV_p(p)  # Enthalpy of vapor
        self.uf = self.steamTable.uL_p(p)  # Internal energy of liquid
        self.ug = self.steamTable.uV_p(p)  # Internal energy of vapor
        self.sf = self.steamTable.sL_p(p)  # Entropy of liquid
        self.sg = self.steamTable.sV_p(p)  # Entropy of vapor
        # Calculate differences between vapor and liquid properties
        self.vgf = self.vg - self.vf
        self.hgf = self.hg - self.hf
        self.sgf = self.sg - self.sf
        self.ugf = self.ug - self.uf

# Class to represent a thermodynamic state and calculate its properties
class thermoState:
    """
    A class to represent a thermodynamic state and calculate its properties.

    Attributes:
        steamTable (XSteam): Instance of XSteam for thermodynamic calculations.
        region (str): Thermodynamic region (e.g., "two-phase", "super-heated vapor").
        p (float): Pressure (bar).
        t (float): Temperature (C).
        v (float): Specific volume (m^3/kg).
        u (float): Internal energy (kJ/kg).
        h (float): Enthalpy (kJ/kg).
        s (float): Entropy (kJ/kg·K).
        x (float): Quality (fraction, 0 to 1).
    """
    def __init__(self, p=None, t=None, v=None, u=None, h=None, s=None, x=None):
        """
        Initialize the thermoState class.

            Args:
                p (float, optional): Pressure in bar. Defaults to None.
                t (float, optional): Temperature in Celsius. Defaults to None.
                v (float, optional): Specific volume in m^3/kg. Defaults to None.
                u (float, optional): Internal energy in kJ/kg. Defaults to None.
                h (float, optional): Enthalpy in kJ/kg. Defaults to None.
                s (float, optional): Entropy in kJ/kg·K. Defaults to None.
                x (float, optional): Quality (0 to 1). Defaults to None.
        """
        # Initialize steam table with SI units (MKS system)
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)
        self.region = "saturated"
        # Initialize state properties
        self.p = p
        self.t = t
        self.v = v
        self.u = u
        self.h = h
        self.s = s
        self.x = x

    def computeProperties(self):
        """
        Compute the thermodynamic properties based on the region and known properties.
        Updates the state properties (u, h, s, v, x) based on the region.
        """
        if self.region == "two-phase":
            self.u = self.steamTable.uL_p(self.p) + self.x * (self.steamTable.uV_p(self.p) - self.steamTable.uL_p(self.p))
            self.h = self.steamTable.hL_p(self.p) + self.x * (self.steamTable.hV_p(self.p) - self.steamTable.hL_p(self.p))
            self.s = self.steamTable.sL_p(self.p) + self.x * (self.steamTable.sV_p(self.p) - self.steamTable.sL_p(self.p))
            self.v = self.steamTable.vL_p(self.p) + self.x * (self.steamTable.vV_p(self.p) - self.steamTable.vL_p(self.p))
            print(f"State v: {self.v}")  # Debug print
        else:
            # Calculate properties for sub-cooled liquid or super-heated vapor
            self.u = self.steamTable.u_pt(self.p, self.t)
            self.h = self.steamTable.h_pt(self.p, self.t)
            self.s = self.steamTable.s_pt(self.p, self.t)
            self.v = self.steamTable.v_pt(self.p, self.t)
            print(f"State v: {self.v}")  # Debug print
            self.x = 1.0 if self.region == "super-heated vapor" else 0.0

    def setState(self, stProp1, stProp2, stPropVal1, stPropVal2, SI=True):
        """
        Set the thermodynamic state based on two specified properties.

        Args:
            stProp1 (str): First property (e.g., 'p', 't', 'x').
            stProp2 (str): Second property (e.g., 'p', 't', 'x').
            stPropVal1 (float): Value of the first property.
            stPropVal2 (float): Value of the second property.
            SI (bool, optional): Use SI units if True, English units if False. Defaults to True.
        """
        # Initialize steam table with the appropriate unit system
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS if SI else XSteam.UNIT_SYSTEM_FLS)
        # Convert property names to lowercase for consistency
        SP = [stProp1, stProp2]
        SP[0] = SP[0].lower()
        SP[1] = SP[1].lower()
        f1 = float(stPropVal1)
        f2 = float(stPropVal2)

        # Case 1: One of the properties is pressure (p)
        if SP[0] == 'p' or SP[1] == 'p':
            oFlipped = SP[0] != 'p'
            SP1 = SP[0] if oFlipped else SP[1]
            self.p = f1 if not oFlipped else f2
            tSat = self.steamTable.tsat_p(self.p)
            if SP1 == 't':
                self.t = f2 if not oFlipped else f1
                tSat = round(tSat)
                if self.t < tSat or self.t > tSat:
                    self.region = "sub-cooled liquid" if self.t < tSat else "super-heated vapor"
                else:
                    self.region = "two-phase"
                    self.x = 0.5
            elif SP1 == 'v':
                self.v = f2 if not oFlipped else f1
                vf = round(self.steamTable.vL_p(self.p), 5)
                vg = round(self.steamTable.vV_p(self.p), 3)
                if self.v < vf or self.v > vg:
                    self.region = "sub-cooled liquid" if self.v < vf else "super-heated vapor"
                    dt = 1.0 if self.v > vg else -1.0
                    fn1 = lambda T: self.v - self.steamTable.v_pt(self.p, T)
                    self.t = fsolve(fn1, [tSat + dt])[0]
                else:
                    self.region = "two-phase"
                    self.x = (self.v - vf) / (vg - vf)
                    self.t = tSat
            elif SP1 == 'u':
                self.u = f2 if not oFlipped else f1
                uf = round(self.steamTable.uL_p(self.p), 5)
                ug = round(self.steamTable.uV_p(self.p), 3)
                ugf = ug - uf
                if self.u < uf or self.u > ug:
                    self.region = "sub-cooled liquid" if self.u < uf else "super-heated vapor"
                    dt = 1.0 if self.u > ug else -1.0
                    fn3 = lambda T: self.u - self.steamTable.u_pt(self.p, T)
                    self.t = fsolve(fn3, [tSat + dt])[0]
                else:
                    self.region = "two-phase"
                    self.x = (self.u - uf) / ugf
                    self.t = tSat
            elif SP1 == 'h':
                self.h = f2 if not oFlipped else f1
                hf = self.steamTable.hL_p(self.p)
                hg = self.steamTable.hV_p(self.p)
                hgf = hg - hf
                if self.h < hf or self.h > hg:
                    self.region = "sub-cooled liquid" if self.h < hf else "super-heated vapor"
                    self.t = self.steamTable.t_ph(self.p, self.h)
                else:
                    self.region = "two-phase"
                    self.x = (self.h - hf) / hgf
                    self.t = tSat
            elif SP1 == 's':
                self.s = f2 if not oFlipped else f1
                sf = self.steamTable.sL_p(self.p)
                sg = the.steamTable.sV_p(self.p)
                sgf = sg - sf
                if self.s < sf or self.s > sg:
                    self.region = "sub-cooled liquid" if self.s < sf else "super-heated vapor"
                    self.t = self.steamTable.t_ps(self.p, self.s)
                else:
                    self.region = "two-phase"
                    self.x = (self.s - sf) / sgf
                    self.t = tSat
            elif SP1 == 'x':
                self.region = "two-phase"
                self.x = f2 if not oFlipped else f1
                self.t = tSat

        # Case 2: One of the properties is temperature (t)
        elif SP[0] == 't' or SP[1] == 't':
            oFlipped = SP[0] != 't'
            SP1 = SP[0] if oFlipped else SP[1]
            self.t = f1 if not oFlipped else f2
            pSat = self.steamTable.psat_t(self.t)
            if SP1 == 'v':
                self.v = f2 if not oFlipped else f1
                vf = self.steamTable.vL_p(pSat)
                vg = self.steamTable.vV_p(pSat)
                vgf = vg - vf
                if self.v < vf or self.v > vg:
                    self.region = "sub-cooled liquid" if self.v < vf else "super-heated vapor"
                    dp = -0.1 if self.v > vg else 0.1
                    fn3 = lambda P: [self.v - self.steamTable.v_pt(P, self.t)]
                    self.p = fsolve(fn3, [pSat + dp])[0]
                else:
                    self.region = "two-phase"
                    self.x = (self.v - vf) / vgf
                    self.p = pSat
            elif SP1 == 'u':
                self.u = f2 if not oFlipped else f1
                uf = self.steamTable.uL_p(pSat)
                ug = self.steamTable.uV_p(pSat)
                ugf = ug - uf
                if self.u < uf or self.u > ug:
                    self.region = "sub-cooled liquid" if self.u < uf else "super-heated vapor"
                    dp = 0.1 if self.u > ug else -0.1
                    fn8 = lambda P: self.u - self.steamTable.u_pt(P, self.t)
                    self.p = fsolve(fn8, [pSat + dp])[0]
                else:
                    self.region = "two-phase"
                    self.x = (self.u - uf) / ugf
                    self.p = pSat
            elif SP1 == 'h':
                self.h = f2 if not oFlipped else f1
                hf = self.steamTable.hL_p(pSat)
                hg = self.steamTable.hV_p(pSat)
                hgf = hg - hf
                if self.h < hf or self.h > hg:
                    self.region = "sub-cooled liquid" if self.h < hf else "super-heated vapor"
                    self.p = self.steamTable.p_th(self.t, self.h)
                else:
                    self.region = "two-phase"
                    self.p = pSat
                    self.x = (self.h - hf) / hgf
            elif SP1 == 's':
                self.s = f2 if not oFlipped else f1
                sf = self.steamTable.sL_p(pSat)
                sg = self.steamTable.sV_p(pSat)
                sgf = sg - sf
                if self.s < sf or self.s > sg:
                    self.region = "sub-cooled liquid" if self.s < sf else "super-heated vapor"
                    self.p = self.steamTable.p_ts(self.t, self.s)
                else:
                    self.region = "two-phase"
                    self.p = pSat
                    self.x = (self.s - sf) / sgf
            elif SP1 == 'x':
                self.x = f2 if not oFlipped else f1
                self.region = "two-phase"
                self.p = pSat

        # Case 3: One of the properties is specific volume (v)
        elif SP[0] == 'v' or SP[1] == 'v':
            oFlipped = SP[0] != 'v'
            SP1 = SP[0] if oFlipped else SP[1]
            self.v = f1 if not oFlipped else f2
            if SP1 == 'h':
                self.h = f2 if not oFlipped else f1
                def fn12(P):
                    hf = self.steamTable.hL_p(P)
                    hg = self.steamTable.hV_p(P)
                    hgf = hg - hf
                    vf = self.steamTable.vL_p(P)
                    vg = self.steamTable.vV_p(P)
                    vgf = vg - vf
                    if self.between(self.h, hf, hg):
                        self.x = (self.h - hf) / hgf
                        return self.v - (vf + self.x * vgf)
                    return self.v - self.steamTable.v_ph(P, self.h)
                self.p = fsolve(fn12, [1.0])[0]
                vf = self.steamTable.vL_p(self.p)
                vg = self.steamTable.vV_p(self.p)
                tsat = self.steamTable.tsat_p(self.p)
                if self.v < vf or self.v > vg:
                    self.region = "sub-cooled liquid" if self.v < vf else "super-heated vapor"
                    dt = -1 if self.v < vf else 1
                    findtgivenv = lambda t: self.v - self.steamTable.v_pt(self.p, t)
                    self.t = fsolve(findtgivenv, [tsat + dt])[0]
                else:
                    self.region = "two-phase"
                    self.t = tsat
                    self.x = (self.v - vf) / (vg - vf)
            elif SP1 == 'u':
                self.u = f2 if not oFlipped else f1
                def fn13(PT):
                    p, t = PT
                    uf = self.steamTable.uL_p(p)
                    ug = self.steamTable.uV_p(p)
                    ugf = ug - uf
                    vf = self.steamTable.vL_p(p)
                    vg = self.steamTable.vV_p(p)
                    vgf = vg - vf
                    if self.between(self.u, uf, ug):
                        self.t = self.steamTable.tsat_p(p)
                        self.x = (self.u - uf) / ugf
                        return [self.v - (vf + self.x * vgf), 0]
                    return [self.v - self.steamTable.v_pt(p, t), self.u - self.steamTable.u_pt(p, t)]
                props = fsolve(fn13, [1, 100])
                self.p = props[0]
                self.t = props[1]
                uf = self.steamTable.uL_p(self.p)
                ug = self.steamTable.uV_p(self.p)
                if self.u < uf or self.u > ug:
                    self.region = "sub-cooled" if self.u < uf else "super-heated"
                else:
                    self.region = "two-phase"
                    self.x = (self.u - uf) / (ug - uf)
            elif SP1 == 's':
                self.s = f2 if not oFlipped else f1
                def fn14(PT):
                    p, t = PT
                    sf = self.steamTable.sL_p(p)
                    sg = self.steamTable.sV_p(p)
                    sgf = sg - sf
                    vf = self.steamTable.vL_p(p)
                    vg = self.steamTable.vV_p(p)
                    vgf = vg - vf
                    if self.between(self.s, sf, sg):
                        self.x = (self.s - sf) / sgf
                        return [self.v - (vf + self.x * vgf), 0.0]
                    return [self.v - self.steamTable.v_pt(p, t), self.s - self.steamTable.s_pt(p, t)]
                props = fsolve(fn14, [1, self.steamTable.tsat_p(1)])
                self.p = props[0]
                self.t = props[1]
                sf = self.steamTable.sL_p(self.p)
                sg = self.steamTable.sV_p(self.p)
                if self.s < sf or self.s > sg:
                    self.region = "sub-cooled liquid" if self.s < sf else "super-heated vapor"
                else:
                    self.region = "two-phase"
                    self.x = (self.s - sf) / (sg - sf)
            elif SP1 == 'x':
                self.x = f2 if not oFlipped else f1
                self.x = self.clamp(self.x, 0.0, 1.0)
                self.region = "two-phase"
                def fn15(p):
                    vf = self.steamTable.vL_p(p)
                    vg = self.steamTable.vV_p(p)
                    vgf = vg - vf
                    return self.v - (vf + self.x * vgf)
                props = fsolve(fn15, [1])
                self.p = props[0]
                self.t = self.steamTable.tsat_p(self.p)

        # Case 4: One of the properties is enthalpy (h)
        elif SP[0] == 'h' or SP[1] == 'h':
            oFlipped = SP[0] != 'h'
            SP1 = SP[0] if oFlipped else SP[1]
            self.h = f1 if not oFlipped else f2
            if SP1 == 'u':
                self.u = f2 if not oFlipped else f1
                def fn16(PT):
                    p, t = PT
                    uf = self.steamTable.uL_p(p)
                    ug = self.steamTable.uV_p(p)
                    ugf = ug - uf
                    hf = self.steamTable.hL_p(p)
                    hg = self.steamTable.hV_p(p)
                    hgf = hg - hf
                    if self.between(self.u, uf, ug):
                        self.x = (self.u - uf) / ugf
                        return [self.h - (hf + self.x * hgf), 0.0]
                    return [self.h - self.steamTable.h_pt(p, t), self.u - self.steamTable.u_pt(p, t)]
                props = fsolve(fn16, [1, 100])
                self.p = props[0]
                self.t = props[1]
                uf = self.steamTable.uL_p(self.p)
                ug = self.steamTable.uV_p(self.p)
                if self.u < uf or self.u > ug:
                    self.region = "sub-cooled liquid" if self.u < uf else "super-heated vapor"
                else:
                    self.region = "two-phase"
                    self.x = (self.u - uf) / (ug - uf)
            elif SP1 == 's':
                self.s = f2 if not oFlipped else f1
                def fn17(PT):
                    p, t = PT
                    sf = self.steamTable.sL_p(p)
                    sg = self.steamTable.sV_p(p)
                    sgf = sg - sf
                    hf = self.steamTable.hL_p(p)
                    hg = self.steamTable.hV_p(p)
                    hgf = hg - hf
                    if self.between(self.s, sf, sg):
                        self.x = (self.s - sf) / sgf
                        return [self.h - (hf + self.x * hgf), 0.0]
                    return [self.h - self.steamTable.h_pt(p, t), self.s - self.steamTable.s_pt(p, t)]
                props = fsolve(fn17, [1, 100])
                self.p = props[0]
                self.t = props[1]
                sf = self.steamTable.sL_p(self.p)
                sg = self.steamTable.sV_p(self.p)
                if self.s < sf or self.s > sg:
                    self.region = "sub-cooled liquid" if self.s < sf else "super-heated vapor"
                else:
                    self.region = "two-phase"
                    self.x = (self.s - sf) / (sg - sf)
            elif SP1 == 'x':
                self.x = f2 if not oFlipped else f1
                self.x = self.clamp(self.x, 0.0, 1.0)
                self.region = "two-phase"
                def fn18(p):
                    hf = self.steamTable.hL_p(p)
                    hg = self.steamTable.hV_p(p)
                    hgf = hg - hf
                    return self.h - (hf + self.x * hgf)
                props = fsolve(fn18, [1])
                self.p = props[0]
                self.t = self.steamTable.tsat_p(self.p)

        # Case 5: One of the properties is internal energy (u)
        elif SP[0] == 'u' or SP[1] == 'u':
            oFlipped = SP[0] != 'u'
            SP1 = SP[0] if oFlipped else SP[1]
            self.u = f1 if not oFlipped else f2
            if SP1 == 's':
                self.s = f2 if not oFlipped else f1
                def fn19(PT):
                    p, t = PT
                    sf = self.steamTable.sL_p(p)
                    sg = self.steamTable.sV_p(p)
                    sgf = sg - sf
                    uf = self.steamTable.uL_p(p)
                    ug = self.steamTable.uV_p(p)
                    ugf = ug - uf
                    if self.between(self.s, sf, sg):
                        self.x = (self.s - sf) / sgf
                        return [self.u - (uf + self.x * ugf), 0.0]
                    return [self.u - self.steamTable.u_pt(p, t), self.s - self.steamTable.s_pt(p, t)]
                props = fsolve(fn19, [1, 100])
                self.p = props[0]
                self.t = props[1]
                sf = self.steamTable.sL_p(self.p)
                sg = self.steamTable.sV_p(self.p)
                if self.s < sf or self.s > sg:
                    self.region = "sub-cooled liquid" if self.s < sf else "super-heated vapor"
                else:
                    self.region = "two-phase"
                    self.x = (self.s - sf) / (sg - sf)
            elif SP1 == 'x':
                self.x = f2 if not oFlipped else f1
                self.x = self.clamp(self.x, 0, 1)
                self.region = "two-phase"
                def fn20(p):
                    uf = self.steamTable.uL_p(p)
                    ug = self.steamTable.uV_p(p)
                    ugf = ug - uf
                    return self.u - (uf + self.x * ugf)
                props = fsolve(fn20, [1])
                self.p = props[0]
                self.t = self.steamTable.tsat_p(self.p)

        # Case 6: One of the properties is entropy (s)
        elif SP[0] == 's' or SP[1] == 's':
            oFlipped = SP[0] != 's'
            SP1 = SP[0] if oFlipped else SP[1]
            self.s = f1 if not oFlipped else f2
            if SP1 == 'x':
                self.x = f2 if not oFlipped else f1
                self.x = self.clamp(self.x, 0, 1)
                self.region = "two-phase"
                def fn21(p):
                    sf = self.steamTable.sL_p(p)
                    sg = self.steamTable.sV_p(p)
                    sgf = sg - sf
                    return self.s - (sf + self.x * sgf)
                props = fsolve(fn21, [1])
                self.p = props[0]
                self.t = self.steamTable.tsat_p(self.p)

        # Compute all properties after setting the state
        self.computeProperties()

    def clamp(self, x, low, high):
        """
        Clamp a value between a lower and upper bound.

        Args:
            x (float): Value to clamp.
            low (float): Lower bound.
            high (float): Upper bound.

        Returns:
            float: Clamped value.
        """
        if x < low:
            return low
        if x > high:
            return high
        return x

    def between(self, x, low, high):
        """
        Check if a value is between two bounds (inclusive).

        Args:
            x (float): Value to check.
            low (float): Lower bound.
            high (float): Upper bound.

        Returns:
            bool: True if x is between low and high, False otherwise.
        """
        return x >= low and x <= high

# Main window class for the Thermodynamic State Calculator GUI
class main_window(QWidget, Ui__frm_StateCalculator):
    """
    Main window class for the Thermodynamic State Calculator GUI.

    Attributes:
        steamTable (XSteam): Instance of XSteam for thermodynamic calculations.
        currentUnits (str): Current unit system ('SI' or 'EN').
        p_Units (str): Units for pressure.
        t_Units (str): Units for temperature.
        u_Units (str): Units for internal energy.
        h_Units (str): Units for enthalpy.
        s_Units (str): Units for entropy.
        v_Units (str): Units for specific volume.
        state1 (thermoState): Thermodynamic state 1.
        state2 (thermoState): Thermodynamic state 2.
    """
    def __init__(self):
        """
        Initialize the main window and set up the UI.
        """
        super().__init__()
        self.setupUi(self)
        self.SetupSlotsAndSignals()
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)
        self.currentUnits = 'SI'
        self.setUnits()
        self.show()

    def SetupSlotsAndSignals(self):
        """
        Connect UI signals to their respective slots (event handlers).
        """
        self._rdo_English.clicked.connect(self.setUnits)
        self._rdo_SI.clicked.connect(self.setUnits)
        self._cmb_Property1_State1.currentIndexChanged.connect(self.setUnits)
        self._cmb_Property2_State1.currentIndexChanged.connect(self.setUnits)
        self._cmb_Property1_State2.currentIndexChanged.connect(self.setUnits)
        self._cmb_Property2_State2.currentIndexChanged.connect(self.setUnits)
        self._pb_Calculate.clicked.connect(self.calculateProperties)

    def setUnits(self):
        """
        Set the unit system (SI or English) and update UI elements accordingly.
        """
        SI = self._rdo_SI.isChecked()
        newUnits = 'SI' if SI else 'EN'
        UnitChange = self.currentUnits != newUnits
        self.currentUnits = newUnits

        if SI:
            self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)
            self.p_Units = "bar"
            self.t_Units = "C"
            self.u_Units = "kJ/kg"
            self.h_Units = "kJ/kg"
            self.s_Units = "kJ/kg·K"
            self.v_Units = "m^3/kg"
        else:
            self.steamTable = XSteam(XSteam.UNIT_SYSTEM_FLS)
            self.p_Units = "psi"
            self.t_Units = "F"
            self.u_Units = "btu/lb"
            self.h_Units = "btu/lb"
            self.s_Units = "btu/lb·F"
            self.v_Units = "ft^3/lb"

        self._update_units_and_values(self._cmb_Property1_State1, self._le_Property1_State1, self._lbl_Property1_Units_State1, SI, UnitChange)
        self._update_units_and_values(self._cmb_Property2_State1, self._le_Property2_State1, self._lbl_Property2_Units_State1, SI, UnitChange)
        self._update_units_and_values(self._cmb_Property1_State2, self._le_Property1_State2, self._lbl_Property1_Units_State2, SI, UnitChange)
        self._update_units_and_values(self._cmb_Property2_State2, self._le_Property2_State2, self._lbl_Property2_Units_State2, SI, UnitChange)

    def _update_units_and_values(self, comboBox, lineEdit, unitLabel, SI, UnitChange):
        """
        Update the units and values in the UI when the unit system changes.

        Args:
            comboBox (QComboBox): Combo box for selecting the property.
            lineEdit (QLineEdit): Line edit for entering the property value.
            unitLabel (QLabel): Label for displaying the unit.
            SI (bool): True if SI units, False if English units.
            UnitChange (bool): True if the unit system has changed.
        """
        prop = comboBox.currentText()
        value = float(lineEdit.text()) if lineEdit.text() else 0.0
        if 'Pressure' in prop:
            unitLabel.setText(self.p_Units)
            if UnitChange:
                value = value * UC.psi_to_bar if SI else value * UC.bar_to_psi
        elif 'Temperature' in prop:
            unitLabel.setText(self.t_Units)
            if UnitChange:
                value = UC.F_to_C(value) if SI else UC.C_to_F(value)
        elif 'Energy' in prop:
            unitLabel.setText(self.u_Units)
            if UnitChange:
                value = value * UC.btuperlb_to_kJperkg if SI else value * UC.kJperkg_to_btuperlb
        elif 'Enthalpy' in prop:
            unitLabel.setText(self.h_Units)
            if UnitChange:
                value = value * UC.btuperlb_to_kJperkg if SI else value * UC.kJperkg_to_btuperlb
        elif 'Entropy' in prop:
            unitLabel.setText(self.s_Units)
            if UnitChange:
                value = value * UC.btuperlbF_to_kJperkgC if SI else value * UC.kJperkgC_to_btuperlbF
        elif 'Volume' in prop:
            unitLabel.setText(self.v_Units)
            if UnitChange:
                value = value * UC.ft3perlb_to_m3perkg if SI else value * UC.m3perkg_to_ft3perlb
        elif 'Quality' in prop:
            unitLabel.setText("")
        lineEdit.setText("{:0.3f}".format(value))

    def makeLabel(self, state):
        """
        Create a formatted string of state properties for display.

        Args:
            state (thermoState): The thermodynamic state to display.

        Returns:
            str: Formatted string of state properties.
        """
        stProps = "Region = {:}\n".format(state.region)
        stProps += "Pressure = {:0.3f} ({:})\n".format(state.p, self.p_Units)
        stProps += "Temperature = {:0.3f} ({:})\n".format(state.t, self.t_Units)
        stProps += "Internal Energy = {:0.3f} ({:})\n".format(state.u, self.u_Units)
        stProps += "Enthalpy = {:0.3f} ({:})\n".format(state.h, self.h_Units)
        stProps += "Entropy = {:0.3f} ({:})\n".format(state.s, self.s_Units)
        stProps += "Specific Volume = {:0.6f} ({:})\n".format(state.v, self.v_Units)
        stProps += "Quality = {:0.3f}".format(state.x)
        return stProps

    def makeDeltaLabel(self, state1, state2):
        """
        Create a formatted string of the change in properties between two states.

        Args:
            state1 (thermoState): First thermodynamic state.
            state2 (thermoState): Second thermodynamic state.

        Returns:
            str: Formatted string of property changes.
        """
        stDelta = "Property Change:\n"
        stDelta += "ΔP = {:0.3f} ({:})\n".format(state2.p - state1.p, self.p_Units)
        stDelta += "ΔT = {:0.3f} ({:})\n".format(state2.t - state1.t, self.t_Units)
        stDelta += "Δu = {:0.3f} ({:})\n".format(state2.u - state1.u, self.u_Units)
        stDelta += "Δh = {:0.3f} ({:})\n".format(state2.h - state1.h, self.h_Units)
        stDelta += "Δs = {:0.3f} ({:})\n".format(state2.s - state1.s, self.s_Units)
        stDelta += "Δv = {:0.3f} ({:})\n".format(state2.v - state1.v, self.v_Units)
        return stDelta

    def calculateProperties(self):
        """
        Calculate the thermodynamic properties for State 1 and State 2, and their differences.
        Updates the UI with the calculated properties.
        """
        self.state1 = thermoState()
        self.state2 = thermoState()
        SI = self._rdo_SI.isChecked()

        # State 1
        SP1 = [self._cmb_Property1_State1.currentText()[-2:-1].lower(), self._cmb_Property2_State1.currentText()[-2:-1].lower()]
        if SP1[0] == SP1[1]:
            self._lbl_Warning.setText("Warning: State 1 - You cannot specify the same property twice.")
            return
        f1 = [float(self._le_Property1_State1.text()), float(self._le_Property2_State1.text())]
        self.state1.setState(SP1[0], SP1[1], f1[0], f1[1], SI)
        self._lbl_State1_Props.setText(self.makeLabel(self.state1))

        # State 2
        SP2 = [self._cmb_Property1_State2.currentText()[-2:-1].lower(), self._cmb_Property2_State2.currentText()[-2:-1].lower()]
        if SP2[0] == SP2[1]:
            self._lbl_Warning.setText("Warning: State 2 - You cannot specify the same property twice.")
            return
        f2 = [float(self._le_Property1_State2.text()), float(self._le_Property2_State2.text())]
        self.state2.setState(SP2[0], SP2[1], f2[0], f2[1], SI)
        self._lbl_State2_Props.setText(self.makeLabel(self.state2))

        # State Change
        self._lbl_StateChange.setText(self.makeDeltaLabel(self.state1, self.state2))
        self._lbl_Warning.setText("")

def main():
    """
    Main function to launch the application.
    """
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)
    app.aboutToQuit.connect(app.deleteLater)
    main_win = main_window()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
#endregion