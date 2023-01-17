import numpy as np

def vaccine_efficacy(x, ic50):
    return(x/(x + ic50))

def vaccine_efficacy_four_antibodies(x, ic50):
    ve_a = vaccine_efficacy(x[0], ic50)
    ve_b = vaccine_efficacy(x[1], ic50)
    ve_c = vaccine_efficacy(x[2], ic50)
    ve_d = vaccine_efficacy(x[3], ic50)
    singletons = ve_a + ve_b + ve_c + ve_d
    pairs = ve_a*(ve_b + ve_c + ve_d) + ve_b*(ve_c + ve_d) + ve_c*ve_d
    triples = ve_a*(ve_b*(ve_c + ve_d)+ ve_c*ve_d) + ve_b*ve_c*ve_d
    quadruples = ve_a * ve_b * ve_c * ve_d
    return(singletons - pairs + triples - quadruples)

def vaccine_efficacy_one_antibodies(x, ic50):
    ve_a = vaccine_efficacy(x[0], ic50)
    return(ve_a)


def sqrt_diff(ic50, days, ve_data, c_dframe):
    res = 0
    for data in ve_data:
        ve_estimate = np.zeros(len(days))
        for i in range(len(data)):
            antibody_level = c_dframe.loc[days[i] - 1][1:5]
            ve_estimate[i] = vaccine_efficacy_four_antibodies(antibody_level, ic50)
        res += np.linalg.norm(data-ve_estimate[0:len(data)])
    return(res)
