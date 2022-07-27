import utils
import labeling

import pandas as pd

def improve_sequence(vienna, s, phred):
    improved_v = list(vienna)
    approximate_bends = utils.series_jerk_points(s)

    for i in approximate_bends:

        adjustment = 0

        if peak := utils.close_to_peak(phred.diff(), i):
            adjustment = peak
        elif valley := utils.close_to_valley(phred.diff(), i):
            adjustment = valley

        print(f"{i}: {adjustment}")

        if adjustment == 1:
            pair_char = vienna[i]
            improved_v[i + 1] = pair_char
        elif adjustment == -1:
            pair_char = vienna[i]
            improved_v[i - 1] = pair_char

    return "".join(improved_v)
