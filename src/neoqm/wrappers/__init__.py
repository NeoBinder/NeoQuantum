from neoqm.wrappers.pyscf_wrapper import PyscfWrapper
from neoqm.wrappers.orca_wrapper import OrcaWrapper




def update_config(config):
    config["directory"] = config.get("qm_directory")
    pass

    

def wrapper_from_config(config):
    if config.get("engine") == "pyscf":
        wrapper_cls = PyscfWrapper
    elif config.get("engine") == "orca":
        wrapper_cls = OrcaWrapper
    else:
        raise ValueError("Unknown engine {}".format(config.get("engine")))
    update_config(config)
    return wrapper_cls(**config)
    