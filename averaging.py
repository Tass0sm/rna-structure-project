import numpy as np
import pandas as pd
import plotting
import parsing

###############################################################################
#                               generic analysis                              #
###############################################################################

def safe_get_window(s, a, b):
    window = np.full(b - a, np.nan)

    start = max(a, 0)
    end = min(len(s), b)

    for i in range(start, end):
        window[i - a] = s.iloc[i]

    return window

def get_windows(sequence, points, window_min, window_max):

    windows = []

    for p in points:
        window = safe_get_window(sequence, p + window_min, p + window_max)
        windows.append(window)

    return pd.DataFrame(windows, columns=list(range(window_min, window_max)))

###############################################################################
#                              dr kim's analysis                              #
###############################################################################

# every bound is exclusive

def get_all_loops(s):
    diff = s.diff()
    i = 0

    loops = []

    while i < len(s) - 1:
        i, loop = parsing.parse_loop(diff, i)
        loops.append(loop)
        i = parsing.parse_until_positive(diff, i)

    return loops

def get_top_loops(s):
    diff = s.diff()
    i = 0

    loops = []

    while i < len(s) - 2:
        i, loop = parsing.parse_top_loop(s, i)
        loops.append(loop)

    return loops

# a < loop[0], loop[1] < b

def safe_get_loop_window(s, a, loop, b):
    loop_start = loop[0]
    loop_end = loop[1]
    loop_length = loop_end - loop_start
    window_length = b - a - loop_length + 1
    window = np.full(window_length, np.nan)

    start = max(a, 0)
    end = min(len(s), b)

    # left side
    for i in range(start, loop_start):
        window[i - a] = s.iloc[i]

    # loop
    avg = s.iloc[loop_start:loop_end].mean()
    window[loop_start - a] = avg

    # right side
    for i in range(loop_end, end):
        window[i - loop_length - a + 1] = s.iloc[i]

    return window

def get_windows_from_loops(sequence, loops, left_length, right_length):
    windows = []

    for loop in loops:
        a = loop[0] - left_length
        b = loop[1] + right_length

        window = safe_get_loop_window(sequence, a, loop, b)
        windows.append(window)

    return pd.DataFrame(windows, columns=list(range(-left_length, right_length + 1)))

def get_windows_from_structure(sequence, structure, left_length, right_length):
    loops = get_all_loops(structure)
    return get_windows_from_loops(sequence, loops, left_length, right_length)
