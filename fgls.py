import sys
from power_series import *
from utils import *
import json




class FGL:
    def __init__(self, fgl = None, log_fgl = None, exp_fgl = None, name = None, autosave = True, savepath = None):
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
            (0, 0): 1,
            (0, 1): 0,
            (1, 0): 1,
            (1, 1): alpha
        }
        self.x_plus_n_series = {}

        if fgl is not None:
            self.fgl = fgl
            self.accuracy = get_accuracy(self.fgl, alpha)
        if log_fgl is not None:
            self.log_fgl = log_fgl
            self.accuracy = get_accuracy(self.log_fgl, alpha)

            if exp_fgl is None:
                self.calculate_exp()
            else:
                self.exp_fgl = exp_fgl
                
            if fgl is None:
                self.get_fgl_from_log()
        
        
        
        
    
    def do_autosave(self):
        if self.autosave:
            self.save(self.savepath)

    def calculate_exp(self):
        if self.exp_fgl is not None:
            return
        
        print("Calculating exponential of formal group law")
        self.exp_fgl = monic_power_series_inverse(self.log_fgl)
        print("Finished computing exponential of formal group law:")
        print(self.exp_fgl)

        self.do_autosave()
    
    def get_fgl_from_log(self):
        if self.fgl is not None:
            return
        assert self.log_fgl is not None
        
        if self.exp_fgl is None:
            self.calculate_exp()

        print("Computing formal group law")

        log1 = self.log_fgl
        log2 = self.log_fgl.subs(alpha, beta)

        taskbar = WIPBar(
            int(1/2 * (self.accuracy + 1) * (self.accuracy + 2)),
            self.accuracy + 1, "Computing logarithm cross-terms", "Completed computing cross terms"
        )
        mt = mult_trunc([log1]*self.accuracy)

        
        
        cross_terms = [[x] for x in mt]

        for j in range(1, self.accuracy+1):
            for i in range(self.accuracy+1 - j):
                cross_terms[i].append(double_truncate(expand(cross_terms[i][j-1] * log2), alpha, beta))
                
                taskbar.update()
        coeffs = {}
        ret = 0

        taskbar = WIPBar(
            int(1/2 * (self.accuracy + 1) * (self.accuracy + 2)),
            0, "Evaluating formal group law", "Finished evaluating formal group law"
        )
        for deg in range(self.accuracy+1):
            for i in range(deg+1):
                j = deg - i
                coeffs[(i, j)] = self.exp_fgl.coeff(alpha, deg) * binomial(deg, i)
                ret += expand(coeffs[(i, j)] * cross_terms[i][j])

                taskbar.update()
            #ret += expand(the_exp.coeff(alpha, deg) * (log1 + log2)**deg)
        print("Truncating formal group law")
        self.fgl = double_truncate(ret, alpha, beta, self.accuracy, True)
        print("Finished computing formal group law")
        print(self.fgl)

        self.do_autosave()
    
    def get_n_series(self, n, k = 1):
        if (n, k) in self.n_series:
            return self.n_series[(n, k)]
        else:
            print("Calculating [" + str(n) + "]^" + str(k))
            if n == 0:
                self.n_series[(n, k)] = 0
            elif n == 1:
                self.n_series[(1, k)] = alpha**k
            elif k == 0:
                self.n_series[(n, 0)] = 1
                self.do_autosave()
            elif k > 1:
                self.n_series[(n, k)] = mult_trunc([self.get_n_series(n, k-1), self.get_n_series(n, 1)], alpha, M = self.accuracy)[-1]
                self.do_autosave()
            else:
                self.n_series[(n, 1)] = 0
                for deg in range(self.accuracy + 1):
                    for i in range(deg + 1):
                        j = deg - i

                        self.n_series[(n, 1)] += mult_trunc([self.fgl.coeff(alpha, i).coeff(beta, j) * alpha**i, self.get_n_series(n - 1, j)], alpha, M = self.accuracy)[-1]
                #self.n_series[(n, 1)] = truncate(self.n_series[(n, 1)], M = self.accuracy)
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
            obj["fgl"] = {}
            for deg in range(self.accuracy + 1):
                for i in range(deg + 1):
                    j = deg - i
                    obj["fgl"][str(i) + "," + str(j)] = str(self.fgl.coeff(alpha, i).coeff(beta, j))
        if self.log_fgl is not None:
            obj["log_fgl"] = str(self.log_fgl)
        if self.exp_fgl is not None:
            obj["exp_fgl"] = {}
            for deg in range(self.accuracy + 1):
                obj["exp_fgl"][str(deg)] = str(self.exp_fgl.coeff(alpha, deg))
        if self.accuracy is not None:
            obj["accuracy"] = self.accuracy

        obj["n_series"] = {
            str(i) + "," + str(j): str(self.n_series[(i, j)])
            for (i, j) in self.n_series
        }

        obj["x_plus_n_series"] = {
            str(n): str(self.x_plus_n_series[n])
            for n in self.x_plus_n_series
        }

        with open(path, "w") as file:
            json.dump(obj, file, indent = 4)
    
    def load(path):
        print("Loading formal group law")
        with open(path, "r") as file:
            obj = json.load(file)
            name = None
            fgl = None
            log_fgl = None
            exp_fgl = None

            local_dict = {
                "alpha": alpha,
                "beta": beta
            }

            if "name" in obj:
                name = obj["name"]
                print("Read name")
            if "fgl" in obj:
                fgl = 0

                taskbar = WIPBar(
                    len(obj["fgl"]),
                    start_message = "Reading formal group law polynomial",
                    end_message = "Finished reading formal group law polynomial"
                )
                for ij in obj["fgl"]:
                    i, j = [int(x) for x in ij.split(",")]
                    fgl += parse_expr(obj["fgl"][ij], local_dict = local_dict) * alpha**i * beta**j
                    taskbar.update()
            if "log_fgl" in obj:
                log_fgl = parse_expr(obj["log_fgl"], local_dict = local_dict)
                print("Read logarithm of formal group law")
            if "exp_fgl" in obj:
                exp_fgl = 0

                taskbar = WIPBar(
                    len(obj["exp_fgl"]),
                    start_message = "Reading exponential of formal group law",
                    end_message = "Finished reading exponential of formal group law"
                )
                for i in obj["exp_fgl"]:
                    exp_fgl += parse_expr(obj["exp_fgl"][i], local_dict = local_dict) * alpha**int(i)
                    taskbar.update()
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
        print("Loaded formal group law")
        return ret
