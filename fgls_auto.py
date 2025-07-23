import sys
from utils import *
import json
from FPS import *



alpha = symbols("alpha")
ALPHA = FPS.from_polynomial(alpha, alpha)

beta = symbols("beta")
BETA = FPS.from_polynomial(beta, beta)

class FGL2:
    def __init__(self, fgl: FPS = None, log_fgl: FPS = None, exp_fgl: FPS = None, name = None, autosave = True, savepath = None):
        assert fgl or log_fgl

        self.fgl = None
        self.log_fgl = None
        self.exp_fgl = None
        self.name = name

        self.autosave = autosave
        self.savepath = None
        if savepath is None and name is not None:
            self.savepath = "./" + name + ".json"
        assert not self.autosave or self.savepath is not None

        self.n_series = {
            0: FPS.get_const(alpha, 0),
            1: FPS.from_polynomial(alpha, alpha)
        }
        self.beta_plus_n_series = {}

        if fgl is not None:
            self.fgl = fgl
        assert log_fgl is not None
        if log_fgl is not None:
            self.log_fgl = log_fgl
            self.log_fgl.name = "LOG_FGL"
            self.log_fgl.save_powers = True

            self.log_fgl_copy = self.log_fgl.comp([(alpha, BETA)])
            self.log_fgl_copy.name = "LOG_FGL2"
            self.log_fgl_copy.save_powers = True

            if exp_fgl is None:
                self.calculate_exp()
            else:
                self.exp_fgl = exp_fgl
                
            if fgl is None:
                # TODO:
                self.get_fgl_from_log()
                pass
        
        
        
    
    def do_autosave(self):
        if self.autosave:
            self.save(self.savepath)

    def calculate_exp(self):
        if self.exp_fgl is not None:
            return
        
        print("Calculating exponential of formal group law")
        self.exp_fgl = self.log_fgl.comp_inv()
        print("Finished computing exponential of formal group law:")
        print(self.exp_fgl)
    
    def get_fgl_from_log(self):
        self.calculate_exp()
        
        self.fgl = self.exp_fgl.comp([
            (alpha, self.log_fgl + self.log_fgl_copy)
        ], save_terms = True)
        self.fgl.is_symmetric = True
        
    
    def get_n_series(self, n):
        if n in self.n_series:
            return self.n_series[n]
        else:
            # k == 1
            # TODO: use FGL??
            self.n_series[n] = self.exp_fgl.comp([
                (alpha, self.log_fgl*n)
            ], save_terms = False)
            self.do_autosave()
        return self.n_series[n]
    
    def get_beta_plus_n_series(self, n):
        if n in self.beta_plus_n_series:
            return self.beta_plus_n_series[n]
        else:
            if n == 1:
                return self.fgl
            self.beta_plus_n_series[n] = self.fgl.comp([
                (alpha, self.get_n_series(n)),
                (beta, BETA)
            ], save_terms = True)

            self.do_autosave()
            return self.beta_plus_n_series[n]

    def as_parseable_obj(self):
        obj = {}
        if self.name is not None:
            obj["name"] = self.name
        if self.fgl is not None:
            obj["fgl"] = self.fgl.as_parseable_obj()
        if self.log_fgl is not None:
            obj["log_fgl"] = self.log_fgl.as_parseable_obj()
        if self.exp_fgl is not None:
            obj["exp_fgl"] = self.exp_fgl.as_parseable_obj()

        obj["n_series"] = {
            str(n): self.n_series[n].as_parseable_obj()
            for n in self.n_series
        }

        obj["beta_plus_n_series"] = {
            str(n): self.beta_plus_n_series[n].as_parseable_obj()
            for n in self.beta_plus_n_series
        }

        return obj
    
    def save(self, path = None):
        if path == None:
            assert self.name is not None
            path = "./" + self.name + ".json"

        with open(path, "w") as file:
            json.dump(self.as_parseable_obj(), file, indent = 4)
    
    def from_parseable_obj(self, obj):
        if "name" in obj:
            self.name = obj["name"]
            print("Read name")
        if "fgl" in obj:
            self.fgl.from_parseable_obj(obj["fgl"])
        if "log_fgl" in obj:
            self.log_fgl.from_parseable_obj(obj["log_fgl"])
        if "exp_fgl" in obj:
            self.exp_fgl.from_parseable_obj(obj["exp_fgl"])
        if "n_series" in obj:
            for n in obj["n_series"]:
                if not int(n) in self.n_series:
                    self.get_n_series(int(n))
                self.n_series[int(n)].from_parseable_obj(obj["n_series"][n])
        if "beta_plus_n_series" in obj:
            for n in obj["beta_plus_n_series"]:
                if not int(n) in self.beta_plus_n_series:
                    self.get_beta_plus_n_series(int(n))
                self.beta_plus_n_series[int(n)].from_parseable_obj(obj["beta_plus_n_series"][n])
    def load(self, path):
        with open(path, "r") as file:
            obj = json.load(file)

            self.from_parseable_obj(obj)