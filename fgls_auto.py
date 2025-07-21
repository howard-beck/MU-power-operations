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
            (0, 0): FPS.get_const(alpha, 1),
            (0, 1): FPS.get_const(alpha, 0),
            (1, 0): FPS.get_const(alpha, 1),
            (1, 1): FPS.from_polynomial(alpha, alpha)
        }
        self.x_plus_n_series = {}

        if fgl is not None:
            self.fgl = fgl
        if log_fgl is not None:
            self.log_fgl = log_fgl

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

        log_alpha = self.log_fgl
        log_beta = self.log_fgl.comp([(alpha, BETA)])
        
        self.fgl = self.exp_fgl.comp([
            (alpha, log_alpha + log_beta)
        ], save_powers = True)
        self.fgl.is_symmetric = True
        
    
    def get_n_series(self, n, k = 1):
        if (n, k) in self.n_series:
            return self.n_series[(n, k)]
        else:
            print("Calculating [" + str(n) + "]^" + str(k))
            if n == 0:
                self.n_series[(n, k)] = FPS.get_const(alpha, 0)
            elif n == 1:
                self.n_series[(1, k)] = FPS.from_polynomial(alpha**k, alpha)
            elif k == 0:
                self.n_series[(n, 0)] = FPS.get_const(alpha, 0)
                self.do_autosave()
            elif k > 1:
                self.n_series[(n, k)] = self.get_n_series(n, k-1) * self.get_n_series(n)
                self.do_autosave()
            else:
                # k == 1
                # TODO: use FGL??
                self.n_series[(n, 1)] = self.exp_fgl.comp([
                    (alpha, self.log_fgl*n)
                ], save_powers = False)
                self.do_autosave()
            return self.n_series[(n, k)]
    
    def get_x_plus_n_series(self, n):
        if n in self.x_plus_n_series:
            return self.x_plus_n_series[n]
        else:
            taskbar = WIPBar(
                int(1/2 * (self.accuracy + 1) * (self.accuracy + 2)),
                start_message = "Calculating x + [" + str(n) + "] in formal group law",
                end_message = "Finished calculating x + [" + str(n) + "] in formal group law"
            )
            ret = 0
            for deg in range(self.accuracy + 1):
                for i in range(deg + 1):
                    j = deg - i

                    ret += expand(self.fgl.coeff(alpha, i).coeff(beta, j) * x**i * self.get_n_series(n, j))
                    taskbar.update()
            ret2 = 0
            for i in range(self.accuracy + 1):
                for j in range(self.accuracy + 1):
                    ret2 += x**i * alpha**j * ret.coeff(x, i).coeff(alpha, j)
            self.x_plus_n_series[n] = ret2

            self.do_autosave()
            return ret

    def save(self, path = None):
        if path == None:
            assert self.name is not None
            path = "./" + self.name + ".json"

        obj = {}
        if self.name is not None:
            obj["name"] = self.name
        if self.fgl is not None:
            obj["fgl"] = self.fgl.as_parseable_obj()
        if self.log_fgl is not None:
            obj["log_fgl"] = self.log_fgl.as_parseable_obj()
        if self.exp_fgl is not None:
            obj["exp_fgl"] = self.exp_fgl.as_parseable_obj()

        """
        obj["n_series"] = {
            str(i) + "," + str(j): str(self.n_series[(i, j)])
            for (i, j) in self.n_series
        }

        obj["x_plus_n_series"] = {
            str(n): str(self.x_plus_n_series[n])
            for n in self.x_plus_n_series
        }
        """

        with open(path, "w") as file:
            json.dump(obj, file, indent = 4)
    
    def load(self, path):
        with open(path, "r") as file:
            obj = json.load(file)

            if "name" in obj:
                self.name = obj["name"]
                print("Read name")
            if "fgl" in obj:
                self.fgl.from_parseable_obj(obj["fgl"])
            if "log_fgl" in obj:
                self.log_fgl.from_parseable_obj(obj["log_fgl"])
            if "exp_fgl" in obj:
                self.exp_fgl.from_parseable_obj(obj["exp_fgl"])
        """
        ret = FGL(fgl = fgl, log_fgl = log_fgl, exp_fgl = exp_fgl, name = name)
        if "n_series" in obj:
            taskbar = WIPBar(
                len(obj["n_series"]),
                start_message = "Reading n-series and their powers",
                end_message = "Finished reading n-series and their powers"
            )
            for ij in obj["n_series"]:
                i, j = [int(x) for x in ij.split(",")]

                ret.n_series[(i, j)] = parse_expr(obj["n_series"][ij], local_dict = local_dict)
                taskbar.update()
        if "x_plus_n_series" in obj:
            taskbar = WIPBar(
                len(obj["x_plus_n_series"]),
                start_message = "Reading formal group law additions to n-series",
                end_message = "Finished reading formal group law additions to n-series"
            )
            for n in obj["x_plus_n_series"]:
                ret.x_plus_n_series[int(n)] = parse_expr(obj["x_plus_n_series"][n], local_dict = local_dict)
                taskbar.update()
        """
        print("Loaded formal group law")
        #return ret
