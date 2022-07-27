def my_any(l):
    for i, x in zip(range(-1, 2), l):
        if x:
            return i

    return False

def at_peak(diff_series, i):
    return diff_series.iloc[i - 1] > 0 and diff_series.iloc[i] < 0

def close_to_peak(diff_series, i):
    if at_valley(diff_series, i):
        return False

    window = range(i - 1, i + 2)
    return my_any(map(lambda j: at_peak(diff_series, j), window))

def at_valley(diff_series, i):
    return diff_series.iloc[i - 1] < 0 and diff_series.iloc[i] > 0

def close_to_valley(diff_series, i):
    if at_peak(diff_series, i):
        return False

    window = range(i - 1, i + 2)
    return my_any(map(lambda j: at_valley(diff_series, j), window))

# def close_to_valley_or_peak(diff_series, i):
#     return close_to_valley(diff_series, i) or close_to_peak(diff_series, i)

def series_jerk(s):
    return s.diff().diff().fillna(value=0.0)

def series_jerk_points(s):
    jerk = series_jerk(s)

    return list(map(int, jerk[jerk != 0.0].index))
