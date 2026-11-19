from typing import Dict, Any
from halogtf.parameters.io_params import IOParams
from halogtf.parameters.grid_params import GridParams
from halogtf.parameters.init_params import (
    InitParams,
    NFWParams,
    make_init_params,
)
from halogtf.parameters.num_params import NumericalParams
from halogtf.parameters.scatter_params import ScatterParams
from halogtf.parameters.cosmo_params import CosmoParams
from halogtf.parameters.stop_params import StopParams

def _init_param(param_class, arg):
    if arg is None:
        return param_class()
    elif isinstance(arg, dict):
        return param_class(**arg)
    elif isinstance(arg, param_class):
        return arg
    else:
        raise TypeError(f"Expected {param_class.__name__}, dict, or None")

class Config:
    """
    Central container for all static simulation parameters.

    Attributes
    ----------
    io      : IOParams
    grid    : GridParams
    init    : InitParams subclass (e.g., NFWParams, ABGParams)
    scatter : ScatterParams
    num     : NumericalParams
    stop    : StopParams
    cosmo   : CosmoParams

    You can modify the initial profile like:
        config.init = "abg"
        config.init = ("abg", {"alpha": 3.5})
        config.init = ABGParams(alpha=3.5)
    """

    def __init__(
        self,
        io=None,
        grid=None,
        init=None,
        scatter=None,
        num=None,
        stop=None,
        cosmo=None,
    ):
        self.io         = _init_param(IOParams, io)
        self.grid       = _init_param(GridParams, grid)
        self._init      = None
        self.init       = init  # goes through the setter
        self.scatter    = _init_param(ScatterParams, scatter)
        self.num        = _init_param(NumericalParams, num)
        self.stop       = _init_param(StopParams, stop)
        self.cosmo      = _init_param(CosmoParams, cosmo)

    @classmethod
    def from_dict(cls, meta: Dict[str, Dict[str, Any]]) -> "Config":
        """
        Build a Config from the nested dict produced by halogtf.io.read.import_metadata().

        Reconstructs parameters only; State.from_dir restart is not implemented.

        The `meta` dict is expected to have sections like:
          "_init", "grid", "io", "num", "stop", "scatter", "cosmo"
        with keys possibly prefixed by underscores (e.g., "_Mvir", "_rmax", ...).
        """
        # Helper: strip leading underscores off keys
        def norm(d: Dict[str, Any]) -> Dict[str, Any]:
            return {k.lstrip("_"): v for k, v in d.items()}
        
        # Robustly fetch sections (allow with or without leading underscore)
        def get_section(name: str) -> Dict[str, Any]:
            if name in meta:
                return meta[name]
            alt = "_" + name
            if alt in meta:
                return meta[alt]
            raise KeyError(f"Missing '{name}' section in metadata.")
        
        # Sections
        _init_raw = get_section("init")
        _grid_raw = get_section("grid")
        _io_raw   = get_section("io")
        _num_raw = get_section("num")
        _stop_raw = get_section("stop")
        _scatter_raw  = get_section("scatter")
        _cosmo_raw = get_section("cosmo")

        init_raw = norm(_init_raw)
        grid_raw = norm(_grid_raw)
        io_raw   = norm(_io_raw)
        num_raw = norm(_num_raw)
        stop_raw = norm(_stop_raw)
        scatter_raw  = norm(_scatter_raw)
        cosmo_raw = norm(_cosmo_raw)

        # Build InitParams via existing factory
        profile = init_raw.pop("profile", "nfw")
        init_params = make_init_params(profile, **init_raw)

        # Build the other parameter blocks using constructors
        io_params   = IOParams(**io_raw)
        grid_params = GridParams(**grid_raw)
        num_params = NumericalParams(**num_raw)
        stop_params = StopParams(**stop_raw)
        scatter_params  = ScatterParams(**scatter_raw)
        cosmo_params = CosmoParams(**cosmo_raw)

        return cls(io=io_params, grid=grid_params, init=init_params,
                   scatter=scatter_params, num=num_params, stop=stop_params, cosmo=cosmo_params)

    @property
    def init(self):
        return self._init

    @init.setter
    def init(self, value):
        if value is None:
            self._init = NFWParams()
        elif isinstance(value, InitParams):
            self._init = value
        elif isinstance(value, str):
            self._init = make_init_params(value)
        elif isinstance(value, tuple) and len(value) == 2:
            profile, kwargs = value
            self._init = make_init_params(profile, **kwargs)
        else:
            raise TypeError(
                "init must be None, an InitParams instance, a profile name string, "
                "or a (profile, dict) tuple"
            )

    def __repr__(self):
        return (
            f"Config(\n"
            f"  io={self.io},\n"
            f"  grid={self.grid},\n"
            f"  init={self.init},\n"
            f"  scatter={self.scatter},\n"
            f"  num={self.num},\n"
            f"  stop={self.stop}\n"
            f"  cosmo={self.cosmo}\n"
            f")"
        )
