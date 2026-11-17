import math

class StopParams:
    """
    Conditions under which evolution terminates.

    Attributes
    ----------
    t_halt : float
        Absolute simulation halt time in t_s units. Positive infinity disables it.
    rho0_halt : float
        Innermost-cell density rho[0]/rho_s at which to halt after t > 50 t_s.
        Positive infinity disables this limit.
    """
    def __init__(
            self, 
            t_halt      : float = 1e3,
            rho0_halt   : float = 1500,
    ):
        self._t_halt = None
        self._rho0_halt = None

        self.t_halt = t_halt
        self.rho0_halt = rho0_halt

    @property
    def t_halt(self):
        return self._t_halt

    @t_halt.setter
    def t_halt(self, value):
        if isinstance(value, bool) or math.isnan(value) or value <= 0:
            raise ValueError("t_halt must be positive")
        self._t_halt = float(value)

    @property
    def rho0_halt(self):
        return self._rho0_halt

    @rho0_halt.setter
    def rho0_halt(self, value):
        if isinstance(value, bool) or math.isnan(value) or value <= 0:
            raise ValueError("rho0_halt must be positive")
        self._rho0_halt = float(value)

    def __repr__(self):
        attrs = [
            attr for attr in dir(self)
            if not attr.startswith('_') and not callable(getattr(self, attr))
        ]
        attr_strs = []
        for attr in attrs:
            value = getattr(self, attr)
            attr_strs.append(f"{attr}={repr(value)}")
        return f"{self.__class__.__name__}({', '.join(attr_strs)})"
